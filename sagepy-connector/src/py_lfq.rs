use pyo3::prelude::*;

use std::collections::HashMap;
use sage_core::fdr::picked_precursor;
use sage_core::lfq::{FeatureMap, IntegrationStrategy, LfqSettings, PeakScoringStrategy, PrecursorId, PrecursorRange, build_feature_map, Peak};
use sage_core::lfq::PrecursorId::{Charged, Combined};
use sage_core::scoring::Feature;
use sage_core::spectrum::{MS1Spectra};

use crate::py_database::{PyIndexedDatabase, PyPeptideIx};
use crate::py_retention_alignment::PyAlignment;
use crate::py_scoring::PyFeature;
use crate::py_spectrum::{PyProcessedSpectrum, PyProcessedIMSpectrum};

#[pyclass]
pub struct PyPeak {
    pub inner: Peak,
}

#[pymethods]
impl PyPeak {
    #[new]
    #[pyo3(signature = (rt, spectral_angle, score, q_value, rt_min, rt_max, mobility_min=None, mobility_max=None))]
    pub fn new(
        rt: usize,
        spectral_angle: f64,
        score: f64,
        q_value: f32,
        rt_min: f32,
        rt_max: f32,
        mobility_min: Option<f32>,
        mobility_max: Option<f32>,
    ) -> Self {
        PyPeak {
            inner: Peak {
                rt,
                spectral_angle,
                score,
                q_value,
                rt_min,
                rt_max,
                mobility_min,
                mobility_max,
            },
        }
    }

    #[getter]
    pub fn rt(&self) -> usize {
        self.inner.rt
    }

    #[getter]
    pub fn spectral_angle(&self) -> f64 {
        self.inner.spectral_angle
    }

    #[getter]
    pub fn score(&self) -> f64 {
        self.inner.score
    }

    #[getter]
    pub fn q_value(&self) -> f32 {
        self.inner.q_value
    }
}

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
#[derive(Clone, Eq, PartialEq, Hash)]
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

    #[getter]
    pub fn is_combined(&self) -> bool {
        match self.inner {
            Combined(_) => true,
            _ => false,
        }
    }

    #[getter]
    pub fn is_charged(&self) -> bool {
        match self.inner {
            Charged(_) => true,
            _ => false,
        }
    }

    #[getter]
    pub fn combined(&self) -> Option<PyPeptideIx> {
        match self.inner {
            Combined(id) => Some(PyPeptideIx { inner: id }),
            _ => None,
        }
    }

    #[getter]
    pub fn charged(&self) -> Option<(PyPeptideIx, u8)> {
        match self.inner {
            Charged((id, charge)) => Some((PyPeptideIx { inner: id }, charge)),
            _ => None,
        }
    }

    #[getter]
    pub fn peptide_id(&self) -> PyPeptideIx {
        match self.inner {
            Combined(id) => PyPeptideIx { inner: id },
            Charged((id, _)) => PyPeptideIx { inner: id },
        }
    }

    #[getter]
    pub fn charge(&self) -> Option<u8> {
        match self.inner {
            Charged((_, charge)) => Some(charge),
            _ => None,
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
        mobility_pct_tolerance: f32,
        rt_tolerance: f32,
    ) -> Self {
        PyLfqSettings {
            inner: LfqSettings {
                peak_scoring: peak_scoring.inner.clone(),
                integration: integration.inner.clone(),
                spectral_angle,
                ppm_tolerance,
                rt_tolerance,
                combine_charge_states,
                mobility_pct_tolerance,
            },
        }
    }

    #[getter]
    pub fn peak_scoring_strategy(&self) -> PyPeakScoringStrategy {
        PyPeakScoringStrategy {
            inner: self.inner.peak_scoring.clone(),
        }
    }

    #[getter]
    pub fn integration_strategy(&self) -> PyIntegrationStrategy {
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
        mobility_hi: f32,
        mobility_lo: f32,
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
                mobility_hi,
                mobility_lo,
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

    pub fn get_num_ranges(&self) -> usize {
        self.inner.ranges.len()
    }

    pub fn quantify(
        &self,
        database: &PyIndexedDatabase,
        ms1: Vec<PyProcessedSpectrum>,
        alignments: Vec<PyAlignment>,
    ) -> PyResult<HashMap<(PyPrecursorId, bool), (PyPeak, Vec<f64>)>> {
        let spectra = ms1.into_iter().map(|s| s.inner).collect();
        let alignments_inner = alignments.into_iter().map(|a| a.inner).collect::<Vec<_>>();

        let ms1_enum = MS1Spectra::NoMobility(spectra);
        let mut areas = self.inner.quantify(&database.inner, &ms1_enum, &alignments_inner);

        let _ = picked_precursor(&mut areas);

        let mut result = HashMap::new();
        for ((precursor, is_decoy), (peak, intensities)) in areas {
            let py_precursor = match precursor {
                Combined(id) => PyPrecursorId { inner: Combined(id) },
                Charged((id, z)) => PyPrecursorId { inner: Charged((id, z)) },
            };
            result.insert((py_precursor, is_decoy), (PyPeak { inner: peak }, intensities));
        }

        Ok(result)
    }

    pub fn quantify_with_mobility(
        &self,
        database: &PyIndexedDatabase,
        ms1: Vec<PyProcessedIMSpectrum>,
        alignments: Vec<PyAlignment>,
    ) -> PyResult<HashMap<(PyPrecursorId, bool), (PyPeak, Vec<f64>)>> {
        let spectra = ms1.into_iter().map(|s| s.inner).collect();
        let alignments_inner = alignments.into_iter().map(|a| a.inner).collect::<Vec<_>>();

        let ms1_enum = MS1Spectra::WithMobility(spectra);
        let mut areas = self.inner.quantify(&database.inner, &ms1_enum, &alignments_inner);

        let _ = picked_precursor(&mut areas);

        let mut result = HashMap::new();
        for ((precursor, is_decoy), (peak, intensities)) in areas {
            let py_precursor = match precursor {
                Combined(id) => PyPrecursorId { inner: Combined(id) },
                Charged((id, z)) => PyPrecursorId { inner: Charged((id, z)) },
            };
            result.insert((py_precursor, is_decoy), (PyPeak { inner: peak }, intensities));
        }

        Ok(result)
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

    pub fn get_num_ranges(&self) -> usize {
        self.ranges.len()
    }
}

#[pyfunction]
pub fn py_build_feature_map(
    settings: PyLfqSettings,
    precursor_charge: (u8, u8),
    features: Vec<PyFeature>,
) -> PyFeatureMap {
    let features: Vec<Feature> = features.iter().map(|f| f.inner.clone()).collect();
    let feature_map = build_feature_map(settings.inner, precursor_charge, features.as_slice());
    PyFeatureMap {
        inner: feature_map,
    }
}

#[pyfunction]
pub fn py_build_feature_map_psm(
    settings: PyLfqSettings,
    precursor_charge: (u8, u8),
    features: Vec<PyFeature>,
) -> PyFeatureMap {
    let features: Vec<Feature> = features.iter().map(|f| f.inner.clone()).collect();
    let feature_map = build_feature_map(settings.inner, precursor_charge, features.as_slice());
    PyFeatureMap {
        inner: feature_map,
    }
}

#[pymodule]
pub fn py_lfq(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyPeakScoringStrategy>()?;
    m.add_class::<PyIntegrationStrategy>()?;
    m.add_class::<PyPrecursorId>()?;
    m.add_class::<PyLfqSettings>()?;
    m.add_class::<PyPrecursorRange>()?;
    m.add_class::<PyFeatureMap>()?;
    m.add_class::<PyQuery>()?;
    m.add_class::<PyPeak>()?;
    m.add_function(wrap_pyfunction!(py_build_feature_map, m)?)?;
    m.add_function(wrap_pyfunction!(py_build_feature_map_psm, m)?)?;
    Ok(())
}
