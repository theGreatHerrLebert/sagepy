use std::collections::BTreeMap;
use pyo3::prelude::*;
use qfdrust::dataset::{PeptideSpectrumMatch, PsmDataset};

#[pyclass]
#[derive(Clone)]
pub struct PyPeptideSpectrumMatch {
    pub inner: PeptideSpectrumMatch,
}

#[pymethods]
impl PyPeptideSpectrumMatch {
    #[new]
    fn new(spec_id: String, peptide_id: u32, proteins: Vec<String>, decoy: bool, score: f64, intensity: f64, features: Option<Vec<(String, f64)>>) -> Self {
        PyPeptideSpectrumMatch {
            inner: PeptideSpectrumMatch {
                spec_id,
                peptide_id,
                proteins,
                decoy,
                score,
                intensity,
                features
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
    fn intensity(&self) -> f64 {
        self.inner.intensity
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
        let mut psm_map = std::collections::BTreeMap::new();
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
        self.inner.psm_map.keys().cloned().collect()
    }

    pub fn inverted_index(&self) -> BTreeMap<(u32, bool), Vec<PyPeptideSpectrumMatch>> {
        let mut inverted_index: BTreeMap<(u32, bool), Vec<PyPeptideSpectrumMatch>> = BTreeMap::new();
        for (_, psms) in &self.inner.psm_map {
            for psm in psms {
                inverted_index.entry((psm.peptide_id, psm.decoy)).or_insert(Vec::new()).push(PyPeptideSpectrumMatch { inner: psm.clone() });
            }
        }

        inverted_index
    }
}

#[pymodule]
pub fn qfdr(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyPeptideSpectrumMatch>()?;
    m.add_class::<PyPsmDataset>()?;
    Ok(())
}