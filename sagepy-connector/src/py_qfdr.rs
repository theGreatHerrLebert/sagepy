use std::collections::BTreeMap;
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use qfdrust::dataset::{PeptideSpectrumMatch, PsmDataset, TDCMethod};

#[pyclass]
#[derive(Clone)]
pub struct PyTDCMethod {
    pub inner: TDCMethod,
}

#[pymethods]
impl PyTDCMethod {
    #[new]
    fn new(method: &str) -> Self {
        PyTDCMethod {
            inner: TDCMethod::from_str(method),
        }
    }
    pub fn to_str(&self) -> &str {
        self.inner.to_str()
    }
}

#[pyclass]
#[derive(Clone)]
pub struct PyPeptideSpectrumMatch {
    pub inner: PeptideSpectrumMatch,
}

#[pymethods]
impl PyPeptideSpectrumMatch {
    #[new]
    fn new(spec_id: String, peptide_id: u32, proteins: Vec<String>,
           decoy: bool, score: f64, intensity_ms1: Option<f64>,
           intensity_ms2: Option<f64>, features: Option<Vec<(String, f64)>>,
           q_value: Option<f64>, confidence: Option<f64>) -> Self {
        PyPeptideSpectrumMatch {
            inner: PeptideSpectrumMatch {
                spec_id,
                peptide_id,
                proteins,
                decoy,
                score,
                intensity_ms1,
                intensity_ms2,
                features,
                q_value,
                confidence,
            },
        }
    }

    #[getter]
    fn spec_id(&self) -> String {
        self.inner.spec_id.clone()
    }

    #[getter]
    fn peptide_id(&self) -> u32 {
        self.inner.peptide_id
    }

    #[getter]
    fn proteins(&self) -> Vec<String> {
        self.inner.proteins.clone()
    }

    #[getter]
    fn decoy(&self) -> bool {
        self.inner.decoy
    }

    #[getter]
    fn score(&self) -> f64 {
        self.inner.score
    }

    #[getter]
    fn intensity_ms1(&self) -> Option<f64> {
        self.inner.intensity_ms1
    }

    #[getter]
    fn intensity_ms2(&self) -> Option<f64> {
        self.inner.intensity_ms2
    }

    #[getter]
    fn features(&self) -> Option<Vec<(String, f64)>> {
        self.inner.features.clone()
    }
}

#[pyclass]
pub struct PyPsmDataset {
    pub inner: PsmDataset,
}

#[pymethods]
impl PyPsmDataset {
    #[new]
    fn new(spec_ids: Vec<String>, matches: Vec<Vec<PyPeptideSpectrumMatch>>) -> Self {
        let mut psm_map = BTreeMap::new();
        let inner_matches: Vec<Vec<_>> = matches.into_iter().map(|m| m.iter().map(|m| m.inner.clone()).collect()).collect();
        for (spec_id, matches) in spec_ids.into_iter().zip(inner_matches.into_iter()) {
            psm_map.insert(spec_id, matches);
        }
        PyPsmDataset {
            inner: PsmDataset {
                psm_map,
            },
        }
    }
    pub fn get_spec_psms(&self, spec_id: String) -> PyResult<Vec<PyPeptideSpectrumMatch>> {
        Ok(self.inner.psm_map.get(&spec_id).unwrap().iter().map(|psm| PyPeptideSpectrumMatch { inner: psm.clone() }).collect())
    }

    #[getter]
    pub fn size(&self) -> usize {
        self.inner.size()
    }

    #[getter]
    pub fn keys(&self) -> Vec<String> {
        self.inner.get_spectra_ids()
    }

    pub fn tdc(&self, method: PyTDCMethod) -> Vec<PyPeptideSpectrumMatch> {
        self.inner.tdc(method.inner).iter().map(|psm| PyPeptideSpectrumMatch { inner: psm.clone() }).collect()
    }
}

#[pymodule]
pub fn qfdr(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyPeptideSpectrumMatch>()?;
    m.add_class::<PyPsmDataset>()?;
    m.add_class::<PyTDCMethod>()?;
    Ok(())
}