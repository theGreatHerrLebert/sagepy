use pyo3::prelude::*;
use sage_core::lfq::{FeatureMap, IntegrationStrategy, LfqSettings, PeakScoringStrategy, PrecursorId, PrecursorRange};
use sage_core::lfq::PrecursorId::{Charged, Combined};
use crate::py_database::PyPeptideIx;

#[pyclass]
pub struct PyPeakScoringStrategy {
    pub inner: PeakScoringStrategy,
}
#[pymethods]
impl PyPeakScoringStrategy {
    #[new]
    pub fn new(
        strategy: &str,
    ) -> Self {
        PyPeakScoringStrategy {
            inner: match strategy {
                "retention_time" => PeakScoringStrategy::RetentionTime,
                "spectral_angle" => PeakScoringStrategy::SpectralAngle,
                "intensity" => PeakScoringStrategy::Intensity,
                "hybrid" => PeakScoringStrategy::Hybrid,
                _ => panic!("Invalid peak scoring strategy"),
            },
        }
    }
    #[getter]
    pub fn strategy(&self) -> String {
        match self.inner {
            PeakScoringStrategy::RetentionTime => "retention_time".to_string(),
            PeakScoringStrategy::SpectralAngle => "spectral_angle".to_string(),
            PeakScoringStrategy::Intensity => "intensity".to_string(),
            PeakScoringStrategy::Hybrid => "hybrid".to_string(),
        }
    }
}


#[pyclass]
pub struct PyIntegrationStrategy {
    pub inner: IntegrationStrategy,
}
#[pymethods]
impl PyIntegrationStrategy {
    #[new]
    pub fn new(
        strategy: &str,
    ) -> Self {
        PyIntegrationStrategy {
            inner: match strategy {
                "apex" => IntegrationStrategy::Apex,
                "sum" => IntegrationStrategy::Sum,
                _ => panic!("Invalid integration strategy"),
            },
        }
    }
    #[getter]
    pub fn strategy(&self) -> String {
        match self.inner {
            IntegrationStrategy::Apex => "apex".to_string(),
            IntegrationStrategy::Sum => "sum".to_string(),
        }
    }
}

#[pyclass]
pub struct PyPrecursorId {
    pub inner: PrecursorId,
}
#[pymethods]
impl PyPrecursorId {
    #[new]
    pub fn new(
        id: &PyPeptideIx,
    ) -> Self {
        PyPrecursorId {
            inner: Combined(id.inner)
        }
    }

    #[staticmethod]
    pub fn from_charged(id: &PyPeptideIx, charge: u8) -> Self {
        PyPrecursorId {
            inner: Charged((id.inner, charge)),
        }
    }
}

#[pyclass]
#[derive(Clone)]
pub struct PyLfqSettings {
    pub inner: LfqSettings,
}

#[pymethods]
impl PyLfqSettings {
    #[new]
    pub fn new(
        peak_scoring: &PyPeakScoringStrategy,
        integration: &PyIntegrationStrategy,
        spectral_angle: f64,
        ppm_tolerance: f32,
        combine_charge_states: bool,
    ) -> Self {
        PyLfqSettings {
            inner: LfqSettings {
                peak_scoring: peak_scoring.inner.clone(),
                integration: integration.inner.clone(),
                spectral_angle,
                ppm_tolerance,
                combine_charge_states,
            },
        }
    }

    #[getter]
    pub fn peak_scoring(&self) -> PyPeakScoringStrategy {
        PyPeakScoringStrategy {
            inner: self.inner.peak_scoring.clone(),
        }
    }

    #[getter]
    pub fn integration(&self) -> PyIntegrationStrategy {
        PyIntegrationStrategy {
            inner: self.inner.integration.clone(),
        }
    }

    #[getter]
    pub fn spectral_angle(&self) -> f64 {
        self.inner.spectral_angle
    }

    #[getter]
    pub fn ppm_tolerance(&self) -> f32 {
        self.inner.ppm_tolerance
    }

    #[getter]
    pub fn combine_charge_states(&self) -> bool {
        self.inner.combine_charge_states
    }
}

#[pyclass]
#[derive(Clone)]
pub struct PyPrecursorRange {
    pub inner: PrecursorRange,
}

#[pymethods]
impl PyPrecursorRange {
    #[new]
    pub fn new(
        rt: f32,
        mass_lo: f32,
        mass_hi: f32,
        charge: u8,
        isotope: usize,
        peptide: PyPeptideIx,
        file_id: usize,
        decoy: bool,
    ) -> Self {
        PyPrecursorRange {
            inner: PrecursorRange {
                rt,
                mass_lo,
                mass_hi,
                charge,
                isotope,
                peptide: peptide.inner,
                file_id,
                decoy,
            },
        }
    }

    #[getter]
    pub fn rt(&self) -> f32 {
        self.inner.rt
    }

    #[getter]
    pub fn mass_lo(&self) -> f32 {
        self.inner.mass_lo
    }

    #[getter]
    pub fn mass_hi(&self) -> f32 {
        self.inner.mass_hi
    }

    #[getter]
    pub fn charge(&self) -> u8 {
        self.inner.charge
    }

    #[getter]
    pub fn isotope(&self) -> usize {
        self.inner.isotope
    }

    #[getter]
    pub fn peptide(&self) -> PyPeptideIx {
        PyPeptideIx {
            inner: self.inner.peptide
        }
    }

    #[getter]
    pub fn file_id(&self) -> usize {
        self.inner.file_id
    }

    #[getter]
    pub fn decoy(&self) -> bool {
        self.inner.decoy
    }
}

#[pyclass]
pub struct PyFeatureMap {
    pub inner: FeatureMap,
}

#[pymethods]
impl PyFeatureMap {
    #[new]
    pub fn new(ranges: Vec<PyPrecursorRange>, min_rts: Vec<f32>, bin_size: usize, settings: PyLfqSettings) -> Self {
        PyFeatureMap {
            inner: FeatureMap {
                ranges: ranges.into_iter().map(|r| r.inner).collect(),
                min_rts,
                bin_size,
                settings: settings.inner,
            }
        }
    }

    #[getter]
    pub fn ranges(&self) -> Vec<PyPrecursorRange> {
        self.inner.ranges.iter().map(|r| PyPrecursorRange { inner: r.clone() }).collect()
    }

    #[getter]
    pub fn min_rts(&self) -> Vec<f32> {
        self.inner.min_rts.clone()
    }

    #[getter]
    pub fn bin_size(&self) -> usize {
        self.inner.bin_size
    }

    #[getter]
    pub fn settings(&self) -> PyLfqSettings {
        PyLfqSettings {
            inner: self.inner.settings.clone()
        }
    }
}

#[pyclass]
struct PyQuery {
    ranges: Vec<PyPrecursorRange>,
    page_lo: usize,
    page_hi: usize,
    bin_size: usize,
    min_rt: f32,
    max_rt: f32,
}

#[pymethods]
impl PyQuery {
    #[new]
    pub fn new(ranges: Vec<PyPrecursorRange>, page_lo: usize, page_hi: usize, bin_size: usize, min_rt: f32, max_rt: f32) -> Self {
        PyQuery {
            ranges,
            page_lo,
            page_hi,
            bin_size,
            min_rt,
            max_rt,
        }
    }

    #[getter]
    pub fn ranges(&self) -> Vec<PyPrecursorRange> {
        self.ranges.clone()
    }

    #[getter]
    pub fn page_lo(&self) -> usize {
        self.page_lo
    }

    #[getter]
    pub fn page_hi(&self) -> usize {
        self.page_hi
    }

    #[getter]
    pub fn bin_size(&self) -> usize {
        self.bin_size
    }

    #[getter]
    pub fn min_rt(&self) -> f32 {
        self.min_rt
    }

    #[getter]
    pub fn max_rt(&self) -> f32 {
        self.max_rt
    }
}

#[pymodule]
pub fn lfq(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyPeakScoringStrategy>()?;
    m.add_class::<PyIntegrationStrategy>()?;
    m.add_class::<PyPrecursorId>()?;
    m.add_class::<PyLfqSettings>()?;
    m.add_class::<PyPrecursorRange>()?;
    m.add_class::<PyFeatureMap>()?;
    m.add_class::<PyQuery>()?;
    Ok(())
}
