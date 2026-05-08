use crate::py_scoring::PyPsm;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use qfdrust::psm::Psm;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use sage_core::database::PeptideIx;
use sage_core::ion_series::Kind;
use sage_core::scoring::{Feature, Fragments};
use std::collections::HashMap;
use std::fs::File;
use std::path::Path;

use parquet::file::reader::{FileReader, SerializedFileReader};
use parquet::record::RowAccessor;
use unimod::unimod::{quantized_mass_to_unimod, quanzie_mass};

const PROTON_MASS: f32 = 1.00727646677;

#[derive(Debug, thiserror::Error)]
pub enum SageResultsError {
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
    #[error("parquet error: {0}")]
    Parquet(#[from] parquet::errors::ParquetError),
    #[error("missing column: {0}")]
    MissingColumn(String),
    #[error("invalid row: {0}")]
    InvalidRow(String),
    #[error("thread pool error: {0}")]
    ThreadPool(#[from] rayon::ThreadPoolBuildError),
}

impl From<SageResultsError> for PyErr {
    fn from(value: SageResultsError) -> Self {
        PyValueError::new_err(value.to_string())
    }
}

#[derive(Clone)]
struct ResultsRow {
    peptide_idx: u32,
    psm_id: usize,
    peptide_len: usize,
    scannr: String,
    rank: u32,
    label: i32,
    expmass: f32,
    calcmass: f32,
    charge: u8,
    precursor_ppm: f32,
    isotope_error: f32,
    fragment_ppm: f32,
    hyperscore: f64,
    delta_next: f64,
    delta_best: f64,
    matched_peaks: u32,
    longest_b: u32,
    longest_y: u32,
    longest_y_pct: f32,
    missed_cleavages: u8,
    matched_intensity_pct: f32,
    scored_candidates: u32,
    poisson: f64,
    sage_discriminant_score: f32,
    posterior_error: f32,
    spectrum_q: f32,
    peptide_q: f32,
    protein_q: f32,
    ms2_intensity: f32,
    rt: f32,
    aligned_rt: f32,
    predicted_rt: f32,
    delta_rt_model: f32,
    ion_mobility: f32,
    predicted_mobility: f32,
    delta_mobility: f32,
    proteins: String,
    peptide: String,
}

struct ColumnIndex {
    by_name: HashMap<String, usize>,
}

impl ColumnIndex {
    fn new(reader: &SerializedFileReader<File>) -> Self {
        let by_name = reader
            .metadata()
            .file_metadata()
            .schema_descr()
            .columns()
            .iter()
            .enumerate()
            .map(|(i, c)| (c.name().to_string(), i))
            .collect();
        Self { by_name }
    }

    fn required(&self, name: &str) -> Result<usize, SageResultsError> {
        self.by_name
            .get(name)
            .copied()
            .ok_or_else(|| SageResultsError::MissingColumn(name.to_string()))
    }

    fn optional(&self, name: &str) -> Option<usize> {
        self.by_name.get(name).copied()
    }
}

fn get_string(row: &parquet::record::Row, idx: usize) -> Result<String, SageResultsError> {
    Ok(row.get_string(idx)?.to_string())
}

fn get_f32(row: &parquet::record::Row, idx: usize) -> Result<f32, SageResultsError> {
    Ok(row
        .get_double(idx)
        .or_else(|_| row.get_float(idx).map(|v| v as f64))? as f32)
}

fn get_f64(row: &parquet::record::Row, idx: usize) -> Result<f64, SageResultsError> {
    Ok(row
        .get_double(idx)
        .or_else(|_| row.get_float(idx).map(|v| v as f64))?)
}

fn get_i64(row: &parquet::record::Row, idx: usize) -> Result<i64, SageResultsError> {
    row.get_long(idx)
        .or_else(|_| row.get_int(idx).map(|v| v as i64))
        .or_else(|_| {
            row.get_ulong(idx).and_then(|v| {
                i64::try_from(v).map_err(|_| {
                    parquet::errors::ParquetError::General(format!(
                        "unsigned value {v} does not fit i64"
                    ))
                })
            })
        })
        .map_err(Into::into)
}

fn get_u32(row: &parquet::record::Row, idx: usize, name: &str) -> Result<u32, SageResultsError> {
    let value = get_i64(row, idx)?;
    u32::try_from(value).map_err(|_| {
        SageResultsError::InvalidRow(format!("column {name} value {value} does not fit u32"))
    })
}

fn get_usize(
    row: &parquet::record::Row,
    idx: usize,
    name: &str,
) -> Result<usize, SageResultsError> {
    let value = get_i64(row, idx)?;
    usize::try_from(value).map_err(|_| {
        SageResultsError::InvalidRow(format!("column {name} value {value} does not fit usize"))
    })
}

fn get_u8(row: &parquet::record::Row, idx: usize, name: &str) -> Result<u8, SageResultsError> {
    let value = get_i64(row, idx)?;
    u8::try_from(value).map_err(|_| {
        SageResultsError::InvalidRow(format!("column {name} value {value} does not fit u8"))
    })
}

fn read_results(path: &Path, max_rank: Option<u32>) -> Result<Vec<ResultsRow>, SageResultsError> {
    let file = File::open(path)?;
    let reader = SerializedFileReader::new(file)?;
    let cols = ColumnIndex::new(&reader);

    let i_psm_id = cols.required("psm_id")?;
    let i_peptide_len = cols.required("peptide_len")?;
    let i_scannr = cols.required("scannr")?;
    let i_rank = cols.required("rank")?;
    let i_label = cols.optional("label");
    let i_is_decoy = cols.optional("is_decoy");
    let i_expmass = cols.required("expmass")?;
    let i_calcmass = cols.required("calcmass")?;
    let i_charge = cols.required("charge")?;
    let i_precursor_ppm = cols.required("precursor_ppm")?;
    let i_isotope_error = cols.required("isotope_error")?;
    let i_fragment_ppm = cols.required("fragment_ppm")?;
    let i_hyperscore = cols.required("hyperscore")?;
    let i_delta_next = cols.required("delta_next")?;
    let i_delta_best = cols.required("delta_best")?;
    let i_matched_peaks = cols.required("matched_peaks")?;
    let i_longest_b = cols.required("longest_b")?;
    let i_longest_y = cols.required("longest_y")?;
    let i_longest_y_pct = cols.required("longest_y_pct")?;
    let i_missed_cleavages = cols.required("missed_cleavages")?;
    let i_matched_intensity_pct = cols.required("matched_intensity_pct")?;
    let i_scored_candidates = cols.required("scored_candidates")?;
    let i_poisson = cols.required("poisson")?;
    let i_sage_discriminant_score = cols.required("sage_discriminant_score")?;
    let i_posterior_error = cols.required("posterior_error")?;
    let i_spectrum_q = cols.required("spectrum_q")?;
    let i_peptide_q = cols.required("peptide_q")?;
    let i_protein_q = cols.required("protein_q")?;
    let i_ms2_intensity = cols.required("ms2_intensity")?;
    let i_rt = cols.required("rt")?;
    let i_aligned_rt = cols.required("aligned_rt")?;
    let i_predicted_rt = cols.required("predicted_rt")?;
    let i_delta_rt_model = cols.required("delta_rt_model")?;
    let i_ion_mobility = cols.required("ion_mobility")?;
    let i_predicted_mobility = cols.required("predicted_mobility")?;
    let i_delta_mobility = cols.required("delta_mobility")?;
    let i_proteins = cols.required("proteins")?;
    let i_peptide = cols.required("peptide")?;

    let mut peptide_to_idx: HashMap<String, u32> = HashMap::new();
    let mut rows = Vec::new();

    for row in reader.get_row_iter(None)? {
        let row = row?;
        let rank = get_u32(&row, i_rank, "rank")?;
        if max_rank.is_some_and(|limit| rank > limit) {
            continue;
        }

        let peptide = get_string(&row, i_peptide)?;
        let next_idx = peptide_to_idx.len() as u32;
        let peptide_idx = *peptide_to_idx.entry(peptide.clone()).or_insert(next_idx);

        let label = match (i_is_decoy, i_label) {
            (Some(idx), _) => {
                let decoy = row.get_bool(idx)?;
                if decoy {
                    -1
                } else {
                    1
                }
            }
            (None, Some(idx)) => {
                let label = i32::try_from(get_i64(&row, idx)?).map_err(|_| {
                    SageResultsError::InvalidRow("label does not fit i32".to_string())
                })?;
                label
            }
            (None, None) => {
                return Err(SageResultsError::MissingColumn(
                    "is_decoy or label".to_string(),
                ))
            }
        };

        rows.push(ResultsRow {
            peptide_idx,
            psm_id: get_usize(&row, i_psm_id, "psm_id")?,
            peptide_len: get_usize(&row, i_peptide_len, "peptide_len")?,
            scannr: get_string(&row, i_scannr)?,
            rank,
            label,
            expmass: get_f32(&row, i_expmass)?,
            calcmass: get_f32(&row, i_calcmass)?,
            charge: get_u8(&row, i_charge, "charge")?,
            precursor_ppm: get_f32(&row, i_precursor_ppm)?,
            isotope_error: get_f32(&row, i_isotope_error)?,
            fragment_ppm: get_f32(&row, i_fragment_ppm)?,
            hyperscore: get_f64(&row, i_hyperscore)?,
            delta_next: get_f64(&row, i_delta_next)?,
            delta_best: get_f64(&row, i_delta_best)?,
            matched_peaks: get_u32(&row, i_matched_peaks, "matched_peaks")?,
            longest_b: get_u32(&row, i_longest_b, "longest_b")?,
            longest_y: get_u32(&row, i_longest_y, "longest_y")?,
            longest_y_pct: get_f32(&row, i_longest_y_pct)?,
            missed_cleavages: get_u8(&row, i_missed_cleavages, "missed_cleavages")?,
            matched_intensity_pct: get_f32(&row, i_matched_intensity_pct)?,
            scored_candidates: get_u32(&row, i_scored_candidates, "scored_candidates")?,
            poisson: get_f64(&row, i_poisson)?,
            sage_discriminant_score: get_f32(&row, i_sage_discriminant_score)?,
            posterior_error: get_f32(&row, i_posterior_error)?,
            spectrum_q: get_f32(&row, i_spectrum_q)?,
            peptide_q: get_f32(&row, i_peptide_q)?,
            protein_q: get_f32(&row, i_protein_q)?,
            ms2_intensity: get_f32(&row, i_ms2_intensity)?,
            rt: get_f32(&row, i_rt)?,
            aligned_rt: get_f32(&row, i_aligned_rt)?,
            predicted_rt: get_f32(&row, i_predicted_rt)?,
            delta_rt_model: get_f32(&row, i_delta_rt_model)?,
            ion_mobility: get_f32(&row, i_ion_mobility)?,
            predicted_mobility: get_f32(&row, i_predicted_mobility)?,
            delta_mobility: get_f32(&row, i_delta_mobility)?,
            proteins: get_string(&row, i_proteins)?,
            peptide,
        });
    }

    Ok(rows)
}

fn read_fragments(path: &Path) -> Result<HashMap<usize, Fragments>, SageResultsError> {
    let file = File::open(path)?;
    let reader = SerializedFileReader::new(file)?;
    let cols = ColumnIndex::new(&reader);

    let i_psm_id = cols.required("psm_id")?;
    let i_fragment_type = cols.required("fragment_type")?;
    let i_fragment_charge = cols.required("fragment_charge")?;
    let i_fragment_ordinals = cols.required("fragment_ordinals")?;
    let i_fragment_intensity = cols.required("fragment_intensity")?;
    let i_fragment_mz_calculated = cols.required("fragment_mz_calculated")?;
    let i_fragment_mz_experimental = cols.required("fragment_mz_experimental")?;

    let mut by_psm_id: HashMap<usize, Fragments> = HashMap::new();

    for row in reader.get_row_iter(None)? {
        let row = row?;
        let psm_id = get_usize(&row, i_psm_id, "psm_id")?;
        let fragment_type = get_string(&row, i_fragment_type)?;
        let kind = match fragment_type.as_str() {
            "b" => Kind::B,
            "y" => Kind::Y,
            other => {
                return Err(SageResultsError::InvalidRow(format!(
                    "unsupported fragment_type {other:?}"
                )))
            }
        };
        let fragments = by_psm_id.entry(psm_id).or_default();
        fragments
            .charges
            .push(get_i64(&row, i_fragment_charge)? as i32);
        fragments.kinds.push(kind);
        fragments
            .fragment_ordinals
            .push(get_i64(&row, i_fragment_ordinals)? as i32);
        fragments
            .intensities
            .push(get_f32(&row, i_fragment_intensity)?);
        fragments
            .mz_calculated
            .push(get_f32(&row, i_fragment_mz_calculated)?);
        fragments
            .mz_experimental
            .push(get_f32(&row, i_fragment_mz_experimental)?);
    }

    Ok(by_psm_id)
}

fn mass_to_unimod(mass: f32, specificity: char) -> Option<&'static str> {
    let quantized = quanzie_mass(mass);
    match (specificity, quantized) {
        ('C', q) if q == quanzie_mass(57.0216) => return Some("[UNIMOD:4]"),
        ('M', q) if q == quanzie_mass(15.9949) => return Some("[UNIMOD:35]"),
        ('[', q) if q == quanzie_mass(42.0) => return Some("[UNIMOD:1]"),
        _ => {}
    }
    quantized_mass_to_unimod()
        .get(&quantized)
        .and_then(|candidates| candidates.first().copied())
}

fn parse_sage_peptide(peptide: &str) -> (String, String) {
    let mut unmodified = String::with_capacity(peptide.len());
    let mut modified = String::with_capacity(peptide.len());
    let chars: Vec<char> = peptide.chars().collect();
    let mut i = 0usize;

    while i < chars.len() {
        if chars[i] == '[' && i + 1 < chars.len() && chars[i + 1] == '+' {
            let start = i;
            i += 2;
            let mass_start = i;
            while i < chars.len() && chars[i] != ']' {
                i += 1;
            }
            if i < chars.len() {
                let mass_str: String = chars[mass_start..i].iter().collect();
                let specificity = if start == 0 { '[' } else { chars[start - 1] };
                if let Ok(mass) = mass_str.parse::<f32>() {
                    if let Some(unimod) = mass_to_unimod(mass, specificity) {
                        modified.push_str(unimod);
                    } else {
                        modified.push('[');
                        modified.push('+');
                        modified.push_str(&mass_str);
                        modified.push(']');
                    }
                }
                i += 1;
                if start == 0 && i < chars.len() && chars[i] == '-' {
                    i += 1;
                }
                continue;
            }
            modified.push(chars[start]);
            unmodified.push(chars[start]);
        } else if chars[i] == '-' {
            i += 1;
            continue;
        } else {
            modified.push(chars[i]);
            unmodified.push(chars[i]);
            i += 1;
        }
    }

    (unmodified, modified)
}

fn row_to_psm(
    row: &ResultsRow,
    fragments: Option<Fragments>,
    default_collision_energy: f32,
) -> Psm {
    let feature = Feature {
        peptide_idx: PeptideIx(row.peptide_idx),
        psm_id: row.psm_id,
        peptide_len: row.peptide_len,
        spec_id: row.scannr.clone(),
        file_id: 0,
        rank: row.rank,
        label: row.label,
        expmass: row.expmass,
        calcmass: row.calcmass,
        charge: row.charge,
        rt: row.rt,
        aligned_rt: row.aligned_rt,
        predicted_rt: row.predicted_rt,
        delta_rt_model: row.delta_rt_model,
        ims: row.ion_mobility,
        predicted_ims: row.predicted_mobility,
        delta_ims_model: row.delta_mobility,
        delta_mass: row.precursor_ppm,
        isotope_error: row.isotope_error,
        average_ppm: row.fragment_ppm,
        hyperscore: row.hyperscore,
        delta_next: row.delta_next,
        delta_best: row.delta_best,
        matched_peaks: row.matched_peaks,
        longest_b: row.longest_b,
        longest_y: row.longest_y,
        longest_y_pct: row.longest_y_pct,
        missed_cleavages: row.missed_cleavages,
        matched_intensity_pct: row.matched_intensity_pct,
        scored_candidates: row.scored_candidates,
        poisson: row.poisson,
        discriminant_score: row.sage_discriminant_score,
        posterior_error: row.posterior_error,
        spectrum_q: row.spectrum_q,
        peptide_q: row.peptide_q,
        protein_q: row.protein_q,
        protein_group_q: row.protein_q,
        ms2_intensity: row.ms2_intensity,
        protein_groups: Some(row.proteins.clone()),
        num_protein_groups: row
            .proteins
            .split(';')
            .filter(|protein| !protein.is_empty())
            .count() as u32,
        fragments,
    };

    let proteins = row
        .proteins
        .split(';')
        .filter(|protein| !protein.is_empty())
        .map(ToString::to_string)
        .collect();
    let (sequence, sequence_modified) = parse_sage_peptide(&row.peptide);
    let mut psm = Psm::new(
        row.scannr.clone(),
        row.peptide_idx,
        proteins,
        feature,
        Some(sequence.clone()),
        Some(sequence_modified.clone()),
        Some(sequence),
        Some(sequence_modified),
        None,
        Some(row.ms2_intensity),
        Some(default_collision_energy),
        Some(default_collision_energy),
        None,
        None,
        None,
    );
    psm.mono_mz_calculated =
        Some((row.calcmass + f32::from(row.charge) * PROTON_MASS) / f32::from(row.charge));
    psm
}

pub fn load_sage_psms_from_parquet_impl(
    results_path: &Path,
    matched_fragments_path: Option<&Path>,
    default_collision_energy: f32,
    max_rank: Option<u32>,
    num_threads: usize,
) -> Result<Vec<PyPsm>, SageResultsError> {
    let rows = read_results(results_path, max_rank)?;
    let fragments_by_psm_id = match matched_fragments_path {
        Some(path) => read_fragments(path)?,
        None => HashMap::new(),
    };

    let num_threads = num_threads.max(1);
    let pool = ThreadPoolBuilder::new().num_threads(num_threads).build()?;
    let psms = pool.install(|| {
        rows.par_iter()
            .map(|row| {
                let fragments = fragments_by_psm_id.get(&row.psm_id).cloned();
                PyPsm {
                    inner: row_to_psm(row, fragments, default_collision_energy),
                }
            })
            .collect()
    });
    Ok(psms)
}

#[pyfunction]
#[pyo3(signature = (results_path, matched_fragments_path=None, default_collision_energy=30.0, max_rank=None, num_threads=4))]
pub fn load_sage_psms_from_parquet(
    results_path: String,
    matched_fragments_path: Option<String>,
    default_collision_energy: f32,
    max_rank: Option<u32>,
    num_threads: usize,
) -> PyResult<Vec<PyPsm>> {
    load_sage_psms_from_parquet_impl(
        Path::new(&results_path),
        matched_fragments_path.as_deref().map(Path::new),
        default_collision_energy,
        max_rank,
        num_threads,
    )
    .map_err(Into::into)
}
