use pyo3::prelude::*;
use sage_core::tmt::{Isobaric, Purity, TmtQuant};
use crate::py_scoring::PyFeature;
use crate::py_spectrum::{PyPeak, PyProcessedSpectrum};

#[pyclass]
pub struct PyIsobaric {
    pub inner: Isobaric,
}

#[pymethods]
impl PyIsobaric {
    #[new]
    pub fn new(
        type_name: &str,
    ) -> Self {
        PyIsobaric {
            inner: match type_name {
                "tmt6" => Isobaric::Tmt6,
                "tmt10" => Isobaric::Tmt10,
                "tmt11" => Isobaric::Tmt11,
                "tmt16" => Isobaric::Tmt16,
                "tmt18" => Isobaric::Tmt18,
                _ => panic!("Invalid isobaric type"),
            },
        }
    }
    #[getter]
    pub fn type_name(&self) -> String {
        match self.inner {
            Isobaric::Tmt6 => "tmt6".to_string(),
            Isobaric::Tmt10 => "tmt10".to_string(),
            Isobaric::Tmt11 => "tmt11".to_string(),
            Isobaric::Tmt16 => "tmt16".to_string(),
            Isobaric::Tmt18 => "tmt18".to_string(),
            _ => panic!("Invalid isobaric type"),
        }
    }

    pub fn modification_mass(&self) -> Option<f32> {
        self.inner.modification_mass()
    }
}

#[pyclass]
#[derive(Clone)]
pub struct PyPurity {
    pub inner: Purity,
}

#[pymethods]
impl PyPurity {
    #[new]
    pub fn new(ratio: f32, correct_precursors: usize, incorrect_precursors: usize, ) -> Self {
        PyPurity {
            inner: Purity {
                ratio,
                correct_precursors,
                incorrect_precursors,
            },
        }
    }

    #[getter]
    pub fn ratio(&self) -> f32 {
        self.inner.ratio
    }

    #[getter]
    pub fn correct_precursors(&self) -> usize {
        self.inner.correct_precursors
    }

    #[getter]
    pub fn incorrect_precursors(&self) -> usize {
        self.inner.incorrect_precursors
    }
}

#[pyclass]
#[derive(Clone)]
pub struct PyQuant {
    pub hit: PyFeature,
    pub hit_purity: PyPurity,
    pub spectrum: PyProcessedSpectrum,
    pub chimera: Option<PyFeature>,
    pub chimera_purity: Option<PyPurity>,
    pub intensities: Vec<Option<PyPeak>>,
}

#[pymethods]
impl PyQuant {
    #[new]
    pub fn new(
        hit: PyFeature,
        hit_purity: PyPurity,
        spectrum: PyProcessedSpectrum,
        intensities: Vec<Option<PyPeak>>,
        chimera: Option<PyFeature>,
        chimera_purity: Option<PyPurity>,
    ) -> Self {
        PyQuant {
            hit,
            hit_purity,
            spectrum,
            chimera,
            chimera_purity,
            intensities,
        }
    }

    #[getter]
    pub fn hit(&self) -> PyFeature {
        self.hit.clone()
    }

    #[getter]
    pub fn hit_purity(&self) -> PyPurity {
        self.hit_purity.clone()
    }

    #[getter]
    pub fn spectrum(&self) -> PyProcessedSpectrum {
        self.spectrum.clone()
    }

    #[getter]
    pub fn chimera(&self) -> Option<PyFeature> {
        self.chimera.clone()
    }

    #[getter]
    pub fn chimera_purity(&self) -> Option<PyPurity> {
        self.chimera_purity.clone()
    }

    #[getter]
    pub fn intensities(&self) -> Vec<Option<PyPeak>> {
        self.intensities.clone()
    }
}

#[pyclass]
#[derive(Clone)]
pub struct PyTmtQuant {
    pub inner: TmtQuant,
}

#[pymethods]
impl PyTmtQuant {
    #[new]
    pub fn new(
        spec_id: String,
        file_id: usize,
        ion_injection_time: f32,
        peaks: Vec<f32>
    ) -> Self {
        PyTmtQuant {
            inner: TmtQuant {
                spec_id,
                file_id,
                ion_injection_time,
                peaks,
            },
        }
    }

    #[getter]
    pub fn spec_id(&self) -> String {
        self.inner.spec_id.clone()
    }

    #[getter]
    pub fn file_id(&self) -> usize {
        self.inner.file_id
    }

    #[getter]
    pub fn ion_injection_time(&self) -> f32 {
        self.inner.ion_injection_time
    }

    #[getter]
    pub fn peaks(&self) -> Vec<f32> {
        self.inner.peaks.clone()
    }
}


#[pymodule]
pub fn py_tmt(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyIsobaric>()?;
    m.add_class::<PyPurity>()?;
    m.add_class::<PyQuant>()?;
    m.add_class::<PyTmtQuant>()?;
    Ok(())
}
