use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use sage_core::modification::{validate_mods, InvalidModification, ModificationSpecificity};
use std::collections::HashMap;
use std::str::FromStr;

#[pyclass]
#[derive(Clone, Debug, PartialEq, Hash)]
pub struct PyModificationSpecificity {
    pub inner: ModificationSpecificity,
}

#[pymethods]
impl PyModificationSpecificity {
    #[new]
    pub fn new(s: &str) -> PyResult<Self> {
        match ModificationSpecificity::from_str(s) {
            Ok(m) => Ok(PyModificationSpecificity { inner: m }),
            Err(InvalidModification::Empty) => {
                Err(PyValueError::new_err("Empty modification string"))
            }
            Err(InvalidModification::InvalidResidue(c)) => Err(PyValueError::new_err(format!(
                "Invalid modification string: unrecognized residue ({})",
                c
            ))),
            Err(InvalidModification::TooLong(s)) => Err(PyValueError::new_err(format!(
                "Invalid modification string: {} is too long",
                s
            ))),
        }
    }

    #[getter]
    pub fn as_string(&self) -> String {
        self.inner.to_string()
    }
}

impl Eq for PyModificationSpecificity {}

#[pyfunction]
pub fn py_validate_mods(input: Option<&PyDict>) -> HashMap<PyModificationSpecificity, f32> {
    // unwrap the input
    let input = input.map(|d| d.extract::<HashMap<String, f32>>().unwrap());
    // validate the mods
    let output = validate_mods(input);
    // convert to a py dict
    let py_validated_mods = output
        .iter()
        .map(|(k, v)| (PyModificationSpecificity { inner: k.clone() }, *v))
        .collect::<HashMap<PyModificationSpecificity, f32>>();

    py_validated_mods
}

#[pyfunction]
pub fn py_validate_var_mods(
    input: Option<&PyDict>,
) -> HashMap<PyModificationSpecificity, Vec<f32>> {
    // unwrap the input
    let input = input.map(|d| d.extract::<HashMap<String, Vec<f32>>>().unwrap());
    let mut output: HashMap<PyModificationSpecificity, Vec<f32>> = HashMap::new();

    if let Some(input) = input {
        for (s, mass) in input {
            match ModificationSpecificity::from_str(&s) {
                Ok(m) => {
                    output.insert(PyModificationSpecificity { inner: m }, mass);
                }
                Err(InvalidModification::Empty) => {
                    log::error!("Skipping invalid modification string: empty")
                }
                Err(InvalidModification::InvalidResidue(c)) => {
                    log::error!(
                        "Skipping invalid modification string: unrecognized residue ({})",
                        c
                    )
                }
                Err(InvalidModification::TooLong(s)) => {
                    log::error!("Skipping invalid modification string: {} is too long", s)
                }
            }
        }
    }
    output
}

#[pymodule]
pub fn modification(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyModificationSpecificity>()?;
    m.add_wrapped(wrap_pyfunction!(py_validate_mods))?;
    m.add_wrapped(wrap_pyfunction!(py_validate_var_mods))?;
    Ok(())
}
