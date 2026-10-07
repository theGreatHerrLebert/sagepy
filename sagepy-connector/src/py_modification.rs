use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use sage_core::modification::{
    validate_mods, InvalidModification, ModificationSpecificity, VarModEntry, VariableModification,
};
use std::collections::HashMap;
use std::str::FromStr;

/// A variable modification entry as exchanged with Python: either a bare mass
/// (unrestricted except by `max_variable_mods`) or a `(mass, max_count)` tuple
/// limiting how often this modification may occur on one peptide.
#[derive(Clone, Debug, PartialEq, FromPyObject, IntoPyObject)]
pub enum PyVarModEntry {
    Mass(f32),
    Limited((f32, Option<usize>)),
}

impl From<PyVarModEntry> for VarModEntry {
    fn from(entry: PyVarModEntry) -> Self {
        match entry {
            PyVarModEntry::Mass(mass) | PyVarModEntry::Limited((mass, None)) => VarModEntry::Mass(mass),
            PyVarModEntry::Limited((mass, max_count)) => {
                VarModEntry::Detailed(VariableModification { mass, max_count })
            }
        }
    }
}

impl From<&VarModEntry> for PyVarModEntry {
    fn from(entry: &VarModEntry) -> Self {
        match entry.max_count() {
            None => PyVarModEntry::Mass(entry.mass()),
            Some(max_count) => PyVarModEntry::Limited((entry.mass(), Some(max_count))),
        }
    }
}

#[pyclass(from_py_object)]
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
#[pyo3(signature = (input=None))]
pub fn py_validate_mods(input: Option<&Bound<'_, PyDict>>) -> HashMap<PyModificationSpecificity, f32> {
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
#[pyo3(signature = (input=None))]
pub fn py_validate_var_mods(
    input: Option<&Bound<'_, PyDict>>,
) -> HashMap<PyModificationSpecificity, Vec<PyVarModEntry>> {
    // unwrap the input
    let input = input.map(|d| d.extract::<HashMap<String, Vec<PyVarModEntry>>>().unwrap());
    let mut output: HashMap<PyModificationSpecificity, Vec<PyVarModEntry>> = HashMap::new();

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
pub fn py_modification(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyModificationSpecificity>()?;
    m.add_wrapped(wrap_pyfunction!(py_validate_mods))?;
    m.add_wrapped(wrap_pyfunction!(py_validate_var_mods))?;
    Ok(())
}
