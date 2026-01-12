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
    pub fn sum(compositions: &Bound<'_, PyList>) -> PyResult<PyComposition> {
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
#[derive(Clone, Copy)]
pub struct PyTolerance {
    pub inner: Tolerance,
}

#[pymethods]
impl PyTolerance {
    #[new]
    #[pyo3(signature = (da=None, ppm=None))]
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
        let result = self.inner * rhs as f32;
        Ok(PyTolerance { inner: result })
    }
}

#[pymodule]
pub fn py_mass(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(h2o, m)?)?;
    m.add_function(wrap_pyfunction!(proton, m)?)?;
    m.add_function(wrap_pyfunction!(neutron, m)?)?;
    m.add_function(wrap_pyfunction!(nh3, m)?)?;
    m.add_function(wrap_pyfunction!(py_monoisotopic, m)?)?;
    m.add_class::<PyTolerance>()?;
    m.add_class::<PyComposition>()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tolerance_ppm_bounds() {
        // Test from sage: Tolerance::Ppm(-10.0, 20.0).bounds(1000.0) == (999.99, 1000.02)
        let tol = Tolerance::Ppm(-10.0, 20.0);
        let (lo, hi) = tol.bounds(1000.0);
        assert!((lo - 999.99).abs() < 0.001);
        assert!((hi - 1000.02).abs() < 0.001);
    }

    #[test]
    fn test_tolerance_ppm_bounds_smaller_mass() {
        // Test from sage: Tolerance::Ppm(-10.0, 10.0).bounds(487.0) == (486.99513, 487.00487)
        let tol = Tolerance::Ppm(-10.0, 10.0);
        let (lo, hi) = tol.bounds(487.0);
        assert!((lo - 486.99513).abs() < 0.0001);
        assert!((hi - 487.00487).abs() < 0.0001);
    }

    #[test]
    fn test_tolerance_da_bounds() {
        let tol = Tolerance::Da(-0.5, 0.5);
        let (lo, hi) = tol.bounds(500.0);
        assert!((lo - 499.5).abs() < 0.001);
        assert!((hi - 500.5).abs() < 0.001);
    }

    #[test]
    fn test_tolerance_contains() {
        let tol = Tolerance::Ppm(-10.0, 10.0);
        // 1000.0 +/- 10 ppm = 999.99 to 1000.01
        assert!(tol.contains(1000.0, 1000.005));
        assert!(tol.contains(1000.0, 999.995));
        assert!(!tol.contains(1000.0, 1000.02));
        assert!(!tol.contains(1000.0, 999.97));
    }

    #[test]
    fn test_tolerance_copy() {
        let tol1 = Tolerance::Ppm(-10.0, 10.0);
        let tol2 = tol1; // Copy
        assert_eq!(tol1.bounds(1000.0), tol2.bounds(1000.0));
    }

    #[test]
    fn test_constants() {
        assert!(H2O > 18.0 && H2O < 18.02);
        assert!(PROTON > 1.0 && PROTON < 1.01);
        assert!(NEUTRON > 1.0 && NEUTRON < 1.01);
        assert!(NH3 > 17.0 && NH3 < 17.03);
    }
}
