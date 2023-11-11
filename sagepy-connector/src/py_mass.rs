use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyList;

use sage_core::mass::{
    composition, monoisotopic, Composition, Tolerance, H2O, NEUTRON, NH3, PROTON,
};

#[pyfunction]
fn h2o() -> f32 {
    H2O
}

#[pyfunction]
fn proton() -> f32 {
    PROTON
}

#[pyfunction]
fn neutron() -> f32 {
    NEUTRON
}

#[pyfunction]
fn nh3() -> f32 {
    NH3
}

#[pyfunction]
fn py_monoisotopic(aa: &str) -> PyResult<f32> {
    if aa.len() == 1 && aa.chars().next().unwrap().is_ascii_uppercase() {
        let aa_u8 = aa.as_bytes()[0];
        Ok(monoisotopic(aa_u8))
    } else {
        Err(PyErr::new::<PyValueError, _>(
            "Input must be a single uppercase ASCII character.",
        ))
    }
}

#[pyclass]
#[derive(Clone)]
pub struct PyComposition {
    inner: Composition,
}

#[pymethods]
impl PyComposition {
    #[new]
    pub fn new(carbon: u16, sulfur: u16) -> Self {
        PyComposition {
            inner: Composition::new(carbon, 0, sulfur),
        }
    }

    // Exposing fields for Python access
    #[getter]
    pub fn carbon(&self) -> u16 {
        self.inner.carbon
    }

    #[getter]
    pub fn sulfur(&self) -> u16 {
        self.inner.sulfur
    }

    // Static method to sum compositions
    #[staticmethod]
    pub fn sum(compositions: &PyList) -> PyResult<PyComposition> {
        let mut total_composition = Composition::new(0, 0, 0);

        for comp in compositions.iter() {
            let py_comp: PyComposition = comp.extract()?;
            total_composition.carbon += py_comp.inner.carbon;
            total_composition.sulfur += py_comp.inner.sulfur;
        }

        Ok(PyComposition {
            inner: total_composition,
        })
    }

    #[staticmethod]
    fn py_composition(aa: &str) -> PyResult<PyComposition> {
        // Ensure the string is exactly one character long
        if aa.chars().count() == 1 {
            // Extract the first character
            let aa_char = aa.chars().next().unwrap(); // Safe to use unwrap here as we know it has exactly one character
            Ok(PyComposition {
                inner: composition(aa_char as u8),
            })
        } else {
            // Return an error if the string is not a single character
            Err(PyErr::new::<PyValueError, _>(
                "Expected a single character string",
            ))
        }
    }
}

#[pyclass]
#[derive(Clone)]
pub struct PyTolerance {
    pub inner: Tolerance,
}

#[pymethods]
impl PyTolerance {
    #[new]
    fn new(da: Option<(f32, f32)>, ppm: Option<(f32, f32)>) -> PyResult<Self> {
        let tolerance = match (da, ppm) {
            (Some((lo, hi)), None) => Tolerance::Da(lo, hi),
            (None, Some((lo, hi))) => Tolerance::Ppm(lo, hi),
            _ => {
                return Err(PyValueError::new_err(
                    "Provide either da or ppm values, not both.",
                ))
            }
        };

        Ok(PyTolerance { inner: tolerance })
    }

    #[getter]
    fn da(&self) -> Option<(f32, f32)> {
        match self.inner {
            Tolerance::Da(lo, hi) => Some((lo, hi)),
            _ => None,
        }
    }

    #[getter]
    fn ppm(&self) -> Option<(f32, f32)> {
        match self.inner {
            Tolerance::Ppm(lo, hi) => Some((lo, hi)),
            _ => None,
        }
    }

    fn bounds(&self, center: f32) -> (f32, f32) {
        self.inner.bounds(center)
    }

    fn contains(&self, center: f32, target: f32) -> bool {
        self.inner.contains(center, target)
    }

    #[staticmethod]
    fn ppm_to_delta_mass(center: f32, ppm: f32) -> f32 {
        Tolerance::ppm_to_delta_mass(center, ppm)
    }

    fn __mul__(&self, rhs: f64) -> PyResult<Self> {
        let result = self.inner.clone() * rhs as f32;
        Ok(PyTolerance { inner: result })
    }
}

#[pymodule]
pub fn mass(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(h2o, m)?)?;
    m.add_function(wrap_pyfunction!(proton, m)?)?;
    m.add_function(wrap_pyfunction!(neutron, m)?)?;
    m.add_function(wrap_pyfunction!(nh3, m)?)?;
    m.add_function(wrap_pyfunction!(py_monoisotopic, m)?)?;
    m.add_class::<PyTolerance>()?;
    m.add_class::<PyComposition>()?;
    Ok(())
}
