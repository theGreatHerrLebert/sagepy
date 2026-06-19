use mzdata::io::{MZFileReader, MZReader as MzDataReader};
use mzdata::prelude::{IonMobilityMeasure as _, SpectrumLike as _};
use mzdata::spectrum::{
    MultiLayerSpectrum as MzSpectrum, RefPeakDataLevel as MzRefPeakDataLevel,
    SignalContinuity as MzSignalContinuity,
};
use numpy::{IntoPyArray, PyArray1, PyArrayMethods};
use pyo3::exceptions::PyIOError;
use pyo3::prelude::*;
use rayon::prelude::*;

use crate::py_mass::PyTolerance;
use crate::tdf::{BrukerMS1CentoidingConfig, BrukerProcessingConfig, TdfReader};
use sage_core::mass::Tolerance as SageTolerance;
use sage_core::spectrum::{
    Deisotoped, IMPeak, Peak, Precursor, ProcessedSpectrum, RawSpectrum, Representation,
    SpectrumProcessor,
};

type MzDataSpectrum = MzSpectrum<mzdata::mzpeaks::CentroidPeak, mzdata::mzpeaks::DeconvolutedPeak>;

fn to_representation(signal_continuity: MzSignalContinuity) -> Representation {
    match signal_continuity {
        MzSignalContinuity::Profile => Representation::Profile,
        _ => Representation::Centroid,
    }
}

fn to_tolerance(window: &mzdata::spectrum::IsolationWindow) -> Option<SageTolerance> {
    if window.is_empty() {
        None
    } else {
        Some(SageTolerance::Da(window.lower_bound, window.upper_bound))
    }
}

fn to_precursors(spectrum: &MzDataSpectrum) -> Vec<PyPrecursor> {
    let fallback_ion_mobility = spectrum
        .acquisition()
        .iter()
        .find_map(|scan| scan.ion_mobility())
        .map(|value| value as f32);

    spectrum
        .precursor_iter()
        .filter_map(|precursor| {
            let selected_ion = precursor.ions.first()?;
            let charge = selected_ion.charge.and_then(|value| value.try_into().ok());
            let inverse_ion_mobility = selected_ion
                .ion_mobility()
                .or_else(|| fallback_ion_mobility.map(|value| value as f64))
                .map(|value| value as f32);
            let collision_energy =
                (precursor.activation.energy != 0.0).then_some(precursor.activation.energy);

            Some(PyPrecursor {
                inner: Precursor {
                    mz: selected_ion.mz as f32,
                    intensity: Some(selected_ion.intensity),
                    charge,
                    spectrum_ref: precursor.precursor_id.clone(),
                    isolation_window: to_tolerance(&precursor.isolation_window),
                    inverse_ion_mobility,
                },
                collision_energy,
            })
        })
        .collect()
}

fn extract_arrays(spectrum: &MzDataSpectrum) -> PyResult<(Vec<f32>, Vec<f32>, Option<Vec<f32>>)> {
    match spectrum.peaks() {
        MzRefPeakDataLevel::Missing => Ok((Vec::new(), Vec::new(), None)),
        MzRefPeakDataLevel::RawData(arrays) => {
            let mz = arrays
                .mzs()
                .map_err(|err| PyIOError::new_err(err.to_string()))?
                .iter()
                .map(|value| *value as f32)
                .collect();
            let intensity = arrays
                .intensities()
                .map_err(|err| PyIOError::new_err(err.to_string()))?
                .iter()
                .copied()
                .collect();
            let mobility = arrays
                .ion_mobility()
                .ok()
                .map(|(values, _)| values.iter().map(|value| *value as f32).collect());
            Ok((mz, intensity, mobility))
        }
        MzRefPeakDataLevel::Centroid(peaks) => {
            let mz = peaks.iter().map(|peak| peak.mz as f32).collect();
            let intensity = peaks.iter().map(|peak| peak.intensity).collect();
            Ok((mz, intensity, None))
        }
        MzRefPeakDataLevel::Deconvoluted(peaks) => {
            let mz = peaks.iter().map(|peak| peak.mz() as f32).collect();
            let intensity = peaks.iter().map(|peak| peak.intensity).collect();
            Ok((mz, intensity, None))
        }
    }
}

fn to_raw_spectrum(spectrum: &MzDataSpectrum, file_id: usize) -> PyResult<PyRawSpectrum> {
    let (mz, intensity, mobility) = extract_arrays(spectrum)?;
    let precursors = to_precursors(spectrum);
    let collision_energies = precursors
        .iter()
        .map(|precursor| precursor.collision_energy)
        .collect();
    let ion_injection_time = spectrum
        .acquisition()
        .iter()
        .next()
        .map(|scan| scan.injection_time)
        .unwrap_or_default();

    Ok(PyRawSpectrum {
        inner: RawSpectrum {
            file_id,
            ms_level: spectrum.ms_level(),
            id: spectrum.id().to_string(),
            precursors: precursors
                .into_iter()
                .map(|precursor| precursor.inner)
                .collect(),
            representation: to_representation(spectrum.signal_continuity()),
            scan_start_time: spectrum.start_time() as f32,
            ion_injection_time,
            total_ion_current: spectrum.peaks().tic(),
            mz,
            intensity,
            mobility,
        },
        collision_energies,
    })
}

#[pyfunction]
#[pyo3(signature = (path, file_id=0, ms_level=None, bruker_config=None, requires_ms1=false))]
pub fn read_spectra(
    path: &str,
    file_id: usize,
    ms_level: Option<u8>,
    bruker_config: Option<PyBrukerProcessingConfig>,
    requires_ms1: bool,
) -> PyResult<Vec<PyRawSpectrum>> {
    if is_pmsms_path(path) {
        return read_pmsms_spectra(path, file_id, ms_level);
    }
    if is_d_path(path) {
        return read_tdf_spectra(path, file_id, ms_level, bruker_config, requires_ms1);
    }

    let reader =
        MzDataReader::open_path(path).map_err(|err| PyIOError::new_err(err.to_string()))?;
    let mut spectra = Vec::new();

    for spectrum in reader {
        if ms_level.is_some_and(|level| spectrum.ms_level() != level) {
            continue;
        }
        spectra.push(to_raw_spectrum(&spectrum, file_id)?);
    }

    Ok(spectra)
}

#[pyfunction]
#[pyo3(signature = (path, take_top_n=150, min_deisotope_mz=0.0, deisotope=true, file_id=0, num_threads=4, bruker_config=None, requires_ms1=false))]
pub fn read_processed_spectra(
    path: &str,
    take_top_n: usize,
    min_deisotope_mz: f32,
    deisotope: bool,
    file_id: usize,
    num_threads: usize,
    bruker_config: Option<PyBrukerProcessingConfig>,
    requires_ms1: bool,
) -> PyResult<Vec<PyProcessedSpectrum>> {
    let raw = read_raw_spectra(path, file_id, bruker_config, requires_ms1)?;
    let processor = SpectrumProcessor {
        take_top_n,
        min_deisotope_mz,
        deisotope,
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();

    Ok(pool.install(|| {
        raw.into_par_iter()
            .filter(|spectrum| spectrum.inner.ms_level == 2 && !spectrum.inner.precursors.is_empty())
            .map(|spectrum| PyProcessedSpectrum {
                inner: processor.process(spectrum.inner),
                collision_energies: spectrum.collision_energies,
            })
            .collect()
    }))
}

fn read_raw_spectra(
    path: &str,
    file_id: usize,
    bruker_config: Option<PyBrukerProcessingConfig>,
    requires_ms1: bool,
) -> PyResult<Vec<PyRawSpectrum>> {
    if is_pmsms_path(path) {
        return read_pmsms_spectra(path, file_id, None);
    }
    if is_d_path(path) {
        return read_tdf_spectra(path, file_id, None, bruker_config, requires_ms1);
    }

    let reader =
        MzDataReader::open_path(path).map_err(|err| PyIOError::new_err(err.to_string()))?;
    reader
        .map(|spectrum| to_raw_spectrum(&spectrum, file_id))
        .collect()
}

fn is_pmsms_path(path: &str) -> bool {
    let trimmed = path.trim_end_matches('/');
    trimmed.to_ascii_lowercase().ends_with(".pmsms")
}

fn read_pmsms_spectra(
    path: &str,
    file_id: usize,
    ms_level: Option<u8>,
) -> PyResult<Vec<PyRawSpectrum>> {
    let dir = std::path::Path::new(path);
    let raw = crate::pmsms::parse(dir, file_id)
        .map_err(|err| PyIOError::new_err(err.to_string()))?;
    Ok(raw_to_py_raw_spectra(raw, ms_level))
}

fn is_d_path(path: &str) -> bool {
    let trimmed = path.trim_end_matches('/');
    trimmed.to_ascii_lowercase().ends_with(".d")
}

fn read_tdf_spectra(
    path: &str,
    file_id: usize,
    ms_level: Option<u8>,
    bruker_config: Option<PyBrukerProcessingConfig>,
    requires_ms1: bool,
) -> PyResult<Vec<PyRawSpectrum>> {
    let cfg = bruker_config.map(|c| c.inner).unwrap_or_default();
    let raw = TdfReader
        .parse(path, file_id, cfg, requires_ms1)
        .map_err(|err| PyIOError::new_err(err.to_string()))?;
    Ok(tdf_to_py_raw_spectra(raw, ms_level))
}

fn tdf_to_py_raw_spectra(
    raw: Vec<(RawSpectrum, Vec<Option<f32>>)>,
    ms_level: Option<u8>,
) -> Vec<PyRawSpectrum> {
    raw.into_iter()
        .filter(|(spec, _)| ms_level.map_or(true, |level| spec.ms_level == level))
        .map(|(inner, mut collision_energies)| {
            // Defensive: keep CE Vec length aligned with precursor count.
            if collision_energies.len() != inner.precursors.len() {
                collision_energies = vec![None; inner.precursors.len()];
            }
            PyRawSpectrum {
                inner,
                collision_energies,
            }
        })
        .collect()
}

/// Wrap raw spectra that don't carry per-spectrum collision energy
/// (.pmsms, plain MGF, etc.) — fills the CE vec with None.
fn raw_to_py_raw_spectra(raw: Vec<RawSpectrum>, ms_level: Option<u8>) -> Vec<PyRawSpectrum> {
    raw.into_iter()
        .filter(|spec| ms_level.map_or(true, |level| spec.ms_level == level))
        .map(|inner| {
            let collision_energies = vec![None; inner.precursors.len()];
            PyRawSpectrum {
                inner,
                collision_energies,
            }
        })
        .collect()
}

#[pyclass]
#[derive(Clone)]
pub struct PyBrukerProcessingConfig {
    pub inner: BrukerProcessingConfig,
}

#[pymethods]
impl PyBrukerProcessingConfig {
    #[new]
    #[pyo3(signature = (mz_ppm=5.0, ims_pct=3.0))]
    pub fn new(mz_ppm: f32, ims_pct: f32) -> Self {
        let inner = BrukerProcessingConfig {
            ms2: Default::default(),
            ms1: BrukerMS1CentoidingConfig { mz_ppm, ims_pct },
        };
        Self { inner }
    }

    #[getter]
    pub fn mz_ppm(&self) -> f32 {
        self.inner.ms1.mz_ppm
    }

    #[getter]
    pub fn ims_pct(&self) -> f32 {
        self.inner.ms1.ims_pct
    }

    fn __repr__(&self) -> String {
        format!(
            "BrukerProcessingConfig(mz_ppm={}, ims_pct={})",
            self.inner.ms1.mz_ppm, self.inner.ms1.ims_pct
        )
    }
}

#[pyclass]
#[derive(Clone)]
pub struct PyRepresentation {
    pub inner: Representation,
}

#[pymethods]
impl PyRepresentation {
    #[new]
    pub fn new(representation: String) -> Self {
        match representation.as_str() {
            "centroid" => PyRepresentation {
                inner: Representation::Centroid,
            },
            "profile" => PyRepresentation {
                inner: Representation::Profile,
            },
            _ => PyRepresentation {
                inner: Representation::Centroid,
            },
        }
    }

    #[getter]
    pub fn representation_as_string(&self) -> String {
        match self.inner {
            Representation::Centroid => "CENTROID",
            Representation::Profile => "PROFILE",
        }
        .to_string()
    }
}

#[pyclass]
#[derive(Clone)]
pub struct PyProcessedSpectrum {
    pub inner: ProcessedSpectrum<Peak>,
    pub collision_energies: Vec<Option<f32>>,
}

#[pymethods]
impl PyProcessedSpectrum {
    #[new]
    pub fn new(
        level: u8,
        id: String,
        file_id: usize,
        scan_start_time: f32,
        ion_injection_time: f32,
        precursors: Vec<PyPrecursor>,
        peaks: Vec<PyPeak>,
        total_ion_current: f32,
    ) -> Self {
        let collision_energies = precursors
            .iter()
            .map(|p| p.collision_energy)
            .collect::<Vec<_>>();
        PyProcessedSpectrum {
            inner: ProcessedSpectrum {
                level,
                id,
                file_id,
                scan_start_time,
                ion_injection_time,
                precursors: precursors.into_iter().map(|p| p.inner).collect(),
                peaks: peaks.into_iter().map(|p| p.inner).collect(),
                peak_charges: vec![],
                total_ion_current,
            },
            collision_energies,
        }
    }

    #[getter]
    pub fn level(&self) -> u8 {
        self.inner.level
    }

    #[setter]
    pub fn set_level(&mut self, level: u8) {
        self.inner.level = level;
    }

    #[getter]
    pub fn id(&self) -> String {
        self.inner.id.clone()
    }

    #[setter]
    pub fn set_id(&mut self, id: String) {
        self.inner.id = id;
    }

    #[getter]
    pub fn file_id(&self) -> usize {
        self.inner.file_id
    }

    #[setter]
    pub fn set_file_id(&mut self, file_id: usize) {
        self.inner.file_id = file_id;
    }

    #[getter]
    pub fn scan_start_time(&self) -> f32 {
        self.inner.scan_start_time
    }

    #[getter]
    pub fn ion_injection_time(&self) -> f32 {
        self.inner.ion_injection_time
    }

    #[getter]
    pub fn precursors(&self) -> Vec<PyPrecursor> {
        self.inner
            .precursors
            .iter()
            .zip(self.collision_energies.iter())
            .into_iter()
            .map(|p| PyPrecursor {
                inner: p.0.clone(),
                collision_energy: *p.1,
            })
            .collect()
    }

    #[getter]
    pub fn peaks(&self) -> Vec<PyPeak> {
        self.inner
            .peaks
            .clone()
            .into_iter()
            .map(|p| PyPeak { inner: p })
            .collect()
    }

    pub fn calibrate_mz_ppm(&mut self, ppm: f32) {
        for peak in self.inner.peaks.iter_mut() {
            let ppm_mass = peak.mass / 1e6;
            let ppm_error = ppm_mass * ppm;
            peak.mass -= ppm_error;
        }
    }

    #[getter]
    pub fn collision_energies(&self) -> Vec<Option<f32>> {
        self.collision_energies.clone()
    }

    #[getter]
    pub fn total_ion_current(&self) -> f32 {
        self.inner.total_ion_current
    }

    pub fn extract_ms1_precursor(&self) -> Option<(f32, u8)> {
        self.inner.extract_ms1_precursor()
    }

    pub fn in_isolation_window(&self, mz: f32) -> Option<bool> {
        self.inner.in_isolation_window(mz)
    }
}

#[pyclass]
#[derive(Clone)]
pub struct PyProcessedIMSpectrum {
    pub inner: ProcessedSpectrum<IMPeak>,
    pub collision_energies: Vec<Option<f32>>,
}

#[pymethods]
impl PyProcessedIMSpectrum {
    #[new]
    pub fn new(
        level: u8,
        id: String,
        file_id: usize,
        scan_start_time: f32,
        ion_injection_time: f32,
        precursors: Vec<PyPrecursor>,
        peaks: Vec<PyIMPeak>,
        total_ion_current: f32,
    ) -> Self {
        let collision_energies = precursors
            .iter()
            .map(|p| p.collision_energy)
            .collect::<Vec<_>>();

        PyProcessedIMSpectrum {
            inner: ProcessedSpectrum {
                level,
                id,
                file_id,
                scan_start_time,
                ion_injection_time,
                precursors: precursors.into_iter().map(|p| p.inner).collect(),
                peaks: peaks.into_iter().map(|p| p.inner).collect(),
                peak_charges: vec![],
                total_ion_current,
            },
            collision_energies,
        }
    }

    #[getter]
    pub fn level(&self) -> u8 {
        self.inner.level
    }

    #[setter]
    pub fn set_level(&mut self, level: u8) {
        self.inner.level = level;
    }

    #[getter]
    pub fn id(&self) -> String {
        self.inner.id.clone()
    }

    #[setter]
    pub fn set_id(&mut self, id: String) {
        self.inner.id = id;
    }

    #[getter]
    pub fn file_id(&self) -> usize {
        self.inner.file_id
    }

    #[setter]
    pub fn set_file_id(&mut self, file_id: usize) {
        self.inner.file_id = file_id;
    }

    #[getter]
    pub fn scan_start_time(&self) -> f32 {
        self.inner.scan_start_time
    }

    #[getter]
    pub fn ion_injection_time(&self) -> f32 {
        self.inner.ion_injection_time
    }

    #[getter]
    pub fn precursors(&self) -> Vec<PyPrecursor> {
        self.inner
            .precursors
            .iter()
            .zip(self.collision_energies.iter())
            .map(|(p, ce)| PyPrecursor {
                inner: p.clone(),
                collision_energy: *ce,
            })
            .collect()
    }

    #[getter]
    pub fn peaks(&self) -> Vec<PyIMPeak> {
        self.inner
            .peaks
            .clone()
            .into_iter()
            .map(|p| PyIMPeak { inner: p })
            .collect()
    }

    #[getter]
    pub fn collision_energies(&self) -> Vec<Option<f32>> {
        self.collision_energies.clone()
    }

    #[getter]
    pub fn total_ion_current(&self) -> f32 {
        self.inner.total_ion_current
    }
}

#[pyclass]
#[derive(Clone)]
pub struct PyRawSpectrum {
    pub inner: RawSpectrum,
    pub collision_energies: Vec<Option<f32>>,
}

#[pymethods]
impl PyRawSpectrum {
    #[new]
    #[pyo3(signature = (file_id, ms_level, id, precursors, representation, scan_start_time, ion_injection_time, total_ion_current, mz, intensity, mobility=None))]
    pub fn new(
        file_id: usize,
        ms_level: u8,
        id: String,
        precursors: Vec<PyPrecursor>,
        representation: PyRepresentation,
        scan_start_time: f32,
        ion_injection_time: f32,
        total_ion_current: f32,
        mz: &Bound<'_, PyArray1<f32>>,
        intensity: &Bound<'_, PyArray1<f32>>,
        mobility: Option<Vec<f32>>,
    ) -> Self {
        let mz_vec = unsafe { mz.as_array().to_vec() };
        let intensity_vec = unsafe { intensity.as_array().to_vec() };
        let collision_energies = precursors
            .iter()
            .map(|p| p.collision_energy)
            .collect::<Vec<_>>();

        PyRawSpectrum {
            inner: RawSpectrum {
                file_id,
                ms_level,
                id,
                precursors: precursors.into_iter().map(|p| p.inner).collect(),
                representation: representation.inner,
                scan_start_time,
                ion_injection_time,
                total_ion_current,
                mz: mz_vec,
                intensity: intensity_vec,
                mobility,
            },
            collision_energies,
        }
    }

    #[getter]
    pub fn file_id(&self) -> usize {
        self.inner.file_id
    }

    #[getter]
    pub fn ms_level(&self) -> u8 {
        self.inner.ms_level
    }

    #[getter]
    pub fn id(&self) -> String {
        self.inner.id.clone()
    }

    #[getter]
    pub fn precursors(&self) -> Vec<PyPrecursor> {
        self.inner
            .precursors
            .iter()
            .zip(self.collision_energies.iter())
            .into_iter()
            .map(|p| PyPrecursor {
                inner: p.0.clone(),
                collision_energy: *p.1,
            })
            .collect()
    }

    #[getter]
    pub fn representation(&self) -> PyRepresentation {
        PyRepresentation {
            inner: self.inner.representation,
        }
    }

    #[getter]
    pub fn scan_start_time(&self) -> f32 {
        self.inner.scan_start_time
    }

    #[getter]
    pub fn ion_injection_time(&self) -> f32 {
        self.inner.ion_injection_time
    }

    #[getter]
    pub fn total_ion_current(&self) -> f32 {
        self.inner.total_ion_current
    }

    #[getter]
    pub fn mz(&self, py: Python) -> Py<PyArray1<f32>> {
        self.inner.mz.clone().into_pyarray(py).unbind()
    }

    #[getter]
    pub fn intensity(&self, py: Python) -> Py<PyArray1<f32>> {
        self.inner.intensity.clone().into_pyarray(py).unbind()
    }

    #[getter]
    pub fn mobility(&self, py: Python) -> Option<Py<PyArray1<f32>>> {
        self.inner
            .mobility
            .as_ref()
            .map(|mob| mob.clone().into_pyarray(py).unbind())
    }

    pub fn filter_top_n(&self, n: usize) -> PyRawSpectrum {
        if let Some(mobility_vec) = &self.inner.mobility {
            // Define inline struct for 3D peak
            struct Peak3D {
                mass: f32,
                intensity: f32,
                mobility: f32,
            }

            let mz_iter = self.inner.mz.iter();
            let intensity_iter = self.inner.intensity.iter();

            let mut peaks = mz_iter
                .zip(intensity_iter)
                .zip(mobility_vec.iter())
                .map(|((m, i), mob)| Peak3D {
                    mass: *m,
                    intensity: *i,
                    mobility: *mob,
                })
                .collect::<Vec<_>>();

            peaks.sort_by(|a, b| {
                b.intensity
                    .partial_cmp(&a.intensity)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            let mut top_peaks = peaks.into_iter().take(n).collect::<Vec<_>>();
            top_peaks.sort_by(|a, b| {
                a.mass
                    .partial_cmp(&b.mass)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });

            PyRawSpectrum {
                inner: RawSpectrum {
                    file_id: self.inner.file_id,
                    ms_level: self.inner.ms_level,
                    id: self.inner.id.clone(),
                    precursors: self.inner.precursors.clone(),
                    representation: self.inner.representation,
                    scan_start_time: self.inner.scan_start_time,
                    ion_injection_time: self.inner.ion_injection_time,
                    total_ion_current: self.inner.total_ion_current,
                    mz: top_peaks.iter().map(|p| p.mass).collect(),
                    intensity: top_peaks.iter().map(|p| p.intensity).collect(),
                    mobility: Some(top_peaks.iter().map(|p| p.mobility).collect()),
                },
                collision_energies: self.collision_energies.clone(),
            }
        } else {
            // Define inline struct for 2D peak
            struct Peak {
                mass: f32,
                intensity: f32,
            }

            let mut peaks = self
                .inner
                .mz
                .iter()
                .zip(self.inner.intensity.iter())
                .map(|(m, i)| Peak {
                    mass: *m,
                    intensity: *i,
                })
                .collect::<Vec<_>>();

            peaks.sort_by(|a, b| {
                b.intensity
                    .partial_cmp(&a.intensity)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            let mut top_peaks = peaks.into_iter().take(n).collect::<Vec<_>>();
            top_peaks.sort_by(|a, b| {
                a.mass
                    .partial_cmp(&b.mass)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });

            PyRawSpectrum {
                inner: RawSpectrum {
                    file_id: self.inner.file_id,
                    ms_level: self.inner.ms_level,
                    id: self.inner.id.clone(),
                    precursors: self.inner.precursors.clone(),
                    representation: self.inner.representation,
                    scan_start_time: self.inner.scan_start_time,
                    ion_injection_time: self.inner.ion_injection_time,
                    total_ion_current: self.inner.total_ion_current,
                    mz: top_peaks.iter().map(|p| p.mass).collect(),
                    intensity: top_peaks.iter().map(|p| p.intensity).collect(),
                    mobility: None,
                },
                collision_energies: self.collision_energies.clone(),
            }
        }
    }
}

#[pyclass]
#[derive(PartialEq, Copy, Clone, Default, Debug)]
pub struct PyPeak {
    pub inner: Peak,
}

#[pymethods]
impl PyPeak {
    #[new]
    pub fn new(mass: f32, intensity: f32) -> Self {
        PyPeak {
            inner: Peak { mass, intensity },
        }
    }

    #[getter]
    pub fn mass(&self) -> f32 {
        self.inner.mass
    }

    #[getter]
    pub fn intensity(&self) -> f32 {
        self.inner.intensity
    }
}

#[pyclass]
#[derive(PartialEq, Copy, Clone, Default, Debug)]
pub struct PyIMPeak {
    pub inner: IMPeak,
}

#[pymethods]
impl PyIMPeak {
    #[new]
    pub fn new(mass: f32, intensity: f32, mobility: f32) -> Self {
        PyIMPeak {
            inner: IMPeak {
                mass,
                intensity,
                mobility,
            },
        }
    }

    #[getter]
    pub fn mass(&self) -> f32 {
        self.inner.mass
    }

    #[getter]
    pub fn intensity(&self) -> f32 {
        self.inner.intensity
    }

    #[getter]
    pub fn mobility(&self) -> f32 {
        self.inner.mobility
    }
}

#[pyclass]
#[derive(Clone)]
pub struct PySpectrumProcessor {
    pub inner: SpectrumProcessor,
}

#[pymethods]
impl PySpectrumProcessor {
    #[new]
    pub fn new(take_top_n: usize, min_deisotope_mz: f32, deisotope: bool) -> Self {
        PySpectrumProcessor {
            inner: SpectrumProcessor {
                take_top_n,
                min_deisotope_mz,
                deisotope,
            },
        }
    }

    #[getter]
    pub fn take_top_n(&self) -> usize {
        self.inner.take_top_n
    }

    #[getter]
    pub fn deisotope(&self) -> bool {
        self.inner.deisotope
    }

    pub fn process(&self, spectrum: &PyRawSpectrum) -> PyProcessedSpectrum {
        PyProcessedSpectrum {
            inner: self.inner.process(spectrum.inner.clone()),
            collision_energies: spectrum.collision_energies.clone(),
        }
    }

    pub fn process_collection(
        &self,
        spectra: Vec<PyRawSpectrum>,
        num_threads: usize,
    ) -> Vec<PyProcessedSpectrum> {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .unwrap();
        let processor = self.inner.clone();

        pool.install(|| {
            spectra
                .into_par_iter()
                .map(|spectrum| PyProcessedSpectrum {
                    inner: processor.process(spectrum.inner),
                    collision_energies: spectrum.collision_energies,
                })
                .collect()
        })
    }

    pub fn process_with_mobility(&self, spectrum: &PyRawSpectrum) -> PyProcessedIMSpectrum {
        PyProcessedIMSpectrum {
            inner: self.inner.process_with_mobility(spectrum.inner.clone()),
            collision_energies: spectrum.collision_energies.clone(),
        }
    }
}

#[pyclass]
#[derive(PartialEq, Clone, Debug)]
pub struct PyDeisotoped {
    pub inner: Deisotoped,
}

#[pymethods]
impl PyDeisotoped {
    #[new]
    #[pyo3(signature = (mz, intensity, charge=None, envelope=None))]
    pub fn new(mz: f32, intensity: f32, charge: Option<u8>, envelope: Option<usize>) -> Self {
        PyDeisotoped {
            inner: Deisotoped {
                mz,
                intensity,
                charge,
                envelope,
            },
        }
    }

    #[getter]
    pub fn mz(&self) -> f32 {
        self.inner.mz
    }

    #[getter]
    pub fn intensity(&self) -> f32 {
        self.inner.intensity
    }

    #[getter]
    pub fn charge(&self) -> Option<u8> {
        self.inner.charge
    }

    #[getter]
    pub fn envelope(&self) -> Option<usize> {
        self.inner.envelope
    }
}

#[pyclass]
#[derive(Default, Clone, Debug)]
pub struct PyPrecursor {
    pub inner: Precursor,
    pub collision_energy: Option<f32>,
}

#[pymethods]
impl PyPrecursor {
    #[new]
    #[pyo3(signature = (mz, intensity=None, charge=None, spectrum_ref=None, isolation_window=None, inverse_ion_mobility=None, collision_energy=None))]
    pub fn new(
        mz: f32,
        intensity: Option<f32>,
        charge: Option<u8>,
        spectrum_ref: Option<String>,
        isolation_window: Option<PyTolerance>,
        inverse_ion_mobility: Option<f32>,
        collision_energy: Option<f32>,
    ) -> Self {
        PyPrecursor {
            inner: Precursor {
                mz,
                intensity,
                charge,
                spectrum_ref,
                isolation_window: isolation_window.map(|t| t.inner),
                inverse_ion_mobility,
            },
            collision_energy,
        }
    }

    #[getter]
    pub fn mz(&self) -> f32 {
        self.inner.mz
    }

    pub fn calibrate_mz_ppm(&mut self, ppm: f32) {
        let mz_ppm = self.inner.mz / 1e6;
        let ppm_error = mz_ppm * ppm;
        self.inner.mz -= ppm_error;
    }

    #[getter]
    pub fn intensity(&self) -> Option<f32> {
        self.inner.intensity
    }

    #[getter]
    pub fn charge(&self) -> Option<u8> {
        self.inner.charge
    }

    #[getter]
    pub fn spectrum_ref(&self) -> Option<String> {
        self.inner.spectrum_ref.clone()
    }

    #[getter]
    pub fn isolation_window(&self) -> Option<PyTolerance> {
        self.inner
            .isolation_window
            .clone()
            .map(|t| PyTolerance { inner: t })
    }

    #[getter]
    pub fn inverse_ion_mobility(&self) -> Option<f32> {
        self.inner.inverse_ion_mobility
    }

    #[getter]
    pub fn collision_energy(&self) -> Option<f32> {
        self.collision_energy
    }
}

#[pymodule]
pub fn py_spectrum(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(read_spectra, m)?)?;
    m.add_function(wrap_pyfunction!(read_processed_spectra, m)?)?;
    m.add_class::<PyPeak>()?;
    m.add_class::<PyIMPeak>()?;
    m.add_class::<PyDeisotoped>()?;
    m.add_class::<PyPrecursor>()?;
    m.add_class::<PySpectrumProcessor>()?;
    m.add_class::<PyRepresentation>()?;
    m.add_class::<PyRawSpectrum>()?;
    m.add_class::<PyProcessedSpectrum>()?;
    m.add_class::<PyProcessedIMSpectrum>()?;
    m.add_class::<PyBrukerProcessingConfig>()?;
    Ok(())
}
