//! Reader for the ionmaiden pmsms binary format.
//!
//! Expects a directory with the following layout (typically `sage_input.pmsms/`):
//!   pmsms.mmappet/        — columnar binary: 0.bin=tof(u32), 1.bin=intensity(u32)
//!   tof2mz.mmappet/       — 0.bin: f32 array mapping tof_index → m/z
//!   precursors.parquet    — one row per precursor with columns:
//!                           precursor_idx(u64), mz(f64), rt(f64),
//!                           inv_ion_mobility(f64), charges(i64),
//!                           fragment_spectrum_start(u64), fragment_event_cnt(u64)
//!
//! Ported from MatteoLacki/sage `devel/optimizations` branch.

use std::fs::File;
use std::path::Path;

use memmap2::Mmap;
use parquet::file::reader::{FileReader, SerializedFileReader};
use parquet::record::RowAccessor;
use sage_core::spectrum::{Precursor, RawSpectrum, Representation};

#[derive(Debug, thiserror::Error)]
pub enum PmsmsError {
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
    #[error("parquet error: {0}")]
    Parquet(#[from] parquet::errors::ParquetError),
    #[error("missing column: {0}")]
    MissingColumn(&'static str),
    #[error("tof index {0} out of range (tof2mz len={1})")]
    TofOutOfRange(u32, usize),
}

struct PrecursorRow {
    precursor_idx: u64,
    mz: f32,
    rt_minutes: f32,
    iim: f32,
    /// Decoded charge states (e.g. encoded value 234 → [2, 3, 4]).
    charges: Vec<u8>,
    frag_start: u64,
    frag_count: u64,
}

/// Decode the pipeline's digit-concatenated charge encoding into individual charges.
/// e.g. 23 → [2, 3], 234 → [2, 3, 4], 2 → [2].
fn decode_charges(encoded: i64) -> Vec<u8> {
    if encoded <= 0 {
        return vec![];
    }
    encoded
        .to_string()
        .bytes()
        .filter_map(|b| {
            let d = b.wrapping_sub(b'0');
            if d > 0 && d <= 9 {
                Some(d)
            } else {
                None
            }
        })
        .collect()
}

fn read_precursors(parquet_path: &Path) -> Result<Vec<PrecursorRow>, PmsmsError> {
    let file = File::open(parquet_path)?;
    let reader = SerializedFileReader::new(file)?;

    let schema = reader.metadata().file_metadata().schema_descr();
    let col_idx = |name: &'static str| -> Result<usize, PmsmsError> {
        schema
            .columns()
            .iter()
            .position(|c| c.name() == name)
            .ok_or(PmsmsError::MissingColumn(name))
    };

    let i_mz = col_idx("mz")?;
    let i_rt = col_idx("rt")?;
    let i_iim = col_idx("inv_ion_mobility")?;
    let i_charge = col_idx("charges")?;
    let i_start = col_idx("fragment_spectrum_start")?;
    let i_count = col_idx("fragment_event_cnt")?;
    let i_pidx = col_idx("precursor_idx")?;

    let mut rows: Vec<PrecursorRow> = Vec::new();

    let iter = reader.get_row_iter(None)?;
    for row in iter {
        let row = row?;
        let precursor_idx = row.get_ulong(i_pidx)?;
        let mz = row.get_double(i_mz)? as f32;
        let rt_minutes = (row.get_double(i_rt)? / 60.0) as f32;
        let iim = row.get_double(i_iim)? as f32;
        let charge_raw = row.get_long(i_charge).unwrap_or(0);
        let charges = decode_charges(charge_raw);
        let frag_start = row.get_ulong(i_start)?;
        let frag_count = row.get_ulong(i_count)?;

        if frag_count == 0 {
            continue;
        }

        rows.push(PrecursorRow {
            precursor_idx,
            mz,
            rt_minutes,
            iim,
            charges,
            frag_start,
            frag_count,
        });
    }

    Ok(rows)
}

/// Memory-map a binary file and reinterpret its contents as a typed slice.
///
/// # Safety
/// The file must contain valid data for type `T` (correct alignment, correct
/// element count). For the u32/f32 binary arrays produced by the ionmaiden
/// pipeline this holds: mmap returns page-aligned memory and u32/f32 have
/// 4-byte alignment.
unsafe fn mmap_as_slice<T: Copy>(path: &Path) -> Result<(Mmap, &'static [T]), PmsmsError> {
    let file = File::open(path)?;
    let mmap = Mmap::map(&file)?;
    let n = mmap.len() / std::mem::size_of::<T>();
    let ptr = mmap.as_ptr() as *const T;
    let slice = std::slice::from_raw_parts(ptr, n);
    let slice: &'static [T] = std::mem::transmute(slice);
    Ok((mmap, slice))
}

pub fn parse(dir: &Path, file_id: usize) -> Result<Vec<RawSpectrum>, PmsmsError> {
    let precursor_path = dir.join("precursors.parquet");
    let tof_bin = dir.join("tof2mz.mmappet").join("0.bin");
    let frag_tof_bin = dir.join("pmsms.mmappet").join("0.bin");
    let frag_int_bin = dir.join("pmsms.mmappet").join("1.bin");

    let precursors = read_precursors(&precursor_path)?;

    // SAFETY: see mmap_as_slice — files are written as packed arrays of the
    // given numeric types by the ionmaiden pipeline.
    let (_tof2mz_mmap, tof2mz) = unsafe { mmap_as_slice::<f32>(&tof_bin)? };
    let (_frag_tof_mmap, frag_tof) = unsafe { mmap_as_slice::<u32>(&frag_tof_bin)? };
    let (_frag_int_mmap, frag_int) = unsafe { mmap_as_slice::<u32>(&frag_int_bin)? };

    let mut spectra = Vec::with_capacity(precursors.len());

    for p in &precursors {
        let start = p.frag_start as usize;
        let end = start + p.frag_count as usize;

        let tof_slice = &frag_tof[start..end];
        let int_slice = &frag_int[start..end];

        let mut mz_vec = Vec::with_capacity(tof_slice.len());
        for &tof in tof_slice {
            let idx = tof as usize;
            if idx >= tof2mz.len() {
                return Err(PmsmsError::TofOutOfRange(tof, tof2mz.len()));
            }
            mz_vec.push(tof2mz[idx]);
        }

        let int_vec: Vec<f32> = int_slice.iter().map(|&i| i as f32).collect();
        let total_ion_current: f32 = int_vec.iter().sum();

        // One Precursor per charge state — mirrors how the MGF parser handles
        // CHARGE=234+ (regex extracts digits 2, 3, 4 as separate Precursor entries).
        let precursors: Vec<Precursor> = if p.charges.is_empty() {
            vec![Precursor {
                mz: p.mz,
                intensity: None,
                charge: None,
                spectrum_ref: None,
                isolation_window: None,
                inverse_ion_mobility: Some(p.iim),
            }]
        } else {
            p.charges
                .iter()
                .map(|&c| Precursor {
                    mz: p.mz,
                    intensity: None,
                    charge: Some(c),
                    spectrum_ref: None,
                    isolation_window: None,
                    inverse_ion_mobility: Some(p.iim),
                })
                .collect()
        };

        spectra.push(RawSpectrum {
            file_id,
            ms_level: 2,
            id: format!("precursor_idx={}", p.precursor_idx),
            precursors,
            representation: Representation::Centroid,
            scan_start_time: p.rt_minutes,
            ion_injection_time: 0.0,
            total_ion_current,
            mz: mz_vec,
            intensity: int_vec,
            mobility: None,
        });
    }

    Ok(spectra)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decode_charges_basic() {
        assert_eq!(decode_charges(2), vec![2u8]);
        assert_eq!(decode_charges(23), vec![2u8, 3]);
        assert_eq!(decode_charges(234), vec![2u8, 3, 4]);
        assert_eq!(decode_charges(0), Vec::<u8>::new());
        assert_eq!(decode_charges(-1), Vec::<u8>::new());
    }
}
