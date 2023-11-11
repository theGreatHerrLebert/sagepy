use crate::py_peptide::PyPeptide;
use pyo3::prelude::*;
use sage_core::ion_series::{Ion, Kind};
use sage_core::mass::monoisotopic;

#[pyclass]
#[derive(Clone)]
pub struct PyKind {
    pub inner: Kind,
}

#[pymethods]
impl PyKind {
    #[new]
    fn new(kind: String) -> PyResult<Self> {
        match kind.to_lowercase().as_str() {
            "a" => Ok(PyKind { inner: Kind::A }),
            "b" => Ok(PyKind { inner: Kind::B }),
            "c" => Ok(PyKind { inner: Kind::C }),
            "x" => Ok(PyKind { inner: Kind::X }),
            "y" => Ok(PyKind { inner: Kind::Y }),
            "z" => Ok(PyKind { inner: Kind::Z }),
            _ => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "Invalid Kind value: {}",
                kind
            ))),
        }
    }
    pub fn kind_as_string(&self) -> String {
        format!("{:?}", self.inner)
    }
}

#[pyclass]
pub struct PyIon {
    pub inner: Ion,
}

#[pymethods]
impl PyIon {
    #[new]
    fn new(kind: PyKind, monoisotopic_mass: f32) -> PyResult<Self> {
        let inner_ion = Ion {
            kind: kind.inner, // Conversion from PyKind to Rust Kind
            monoisotopic_mass,
        };
        Ok(PyIon { inner: inner_ion })
    }

    // Getter methods for accessing Ion properties
    #[getter]
    fn kind(&self) -> PyResult<PyKind> {
        Ok(PyKind {
            inner: self.inner.kind,
        })
    }

    #[getter]
    fn monoisotopic_mass(&self) -> PyResult<f32> {
        Ok(self.inner.monoisotopic_mass)
    }
}

#[pyclass]
pub struct PyIonSeries {
    pub kind: PyKind,
    pub cumulative_mass: f32,
    pub peptide: PyPeptide,
}

#[pymethods]
impl PyIonSeries {
    #[new]
    pub fn new(_py: Python, peptide: PyPeptide, kind: PyKind) -> PyResult<Self> {
        const C: f32 = 12.0;
        const O: f32 = 15.994914;
        const H: f32 = 1.007825;
        const PRO: f32 = 1.0072764;
        const N: f32 = 14.003074;
        const NH3: f32 = N + H * 3.0 + PRO;

        let cumulative_mass = match kind.inner {
            Kind::A => peptide.inner.nterm.unwrap_or_default() - (C + O),
            Kind::B => peptide.inner.nterm.unwrap_or_default(),
            Kind::C => peptide.inner.nterm.unwrap_or_default() + NH3,

            Kind::X => {
                peptide.inner.monoisotopic - peptide.inner.nterm.unwrap_or_default()
                    + (C + O - NH3 + N + H)
            }
            Kind::Y => peptide.inner.monoisotopic - peptide.inner.nterm.unwrap_or_default(),
            Kind::Z => peptide.inner.monoisotopic - peptide.inner.nterm.unwrap_or_default() - NH3,
        };

        Ok(Self {
            kind,
            cumulative_mass,
            peptide,
        })
    }

    #[getter]
    fn kind(&self) -> PyResult<PyKind> {
        Ok(self.kind.clone())
    }

    #[getter]
    fn cumulative_mass(&self) -> PyResult<f32> {
        Ok(self.cumulative_mass)
    }

    #[getter]
    fn peptide(&self) -> PyResult<PyPeptide> {
        Ok(self.peptide.clone())
    }

    pub fn get_ion_series(&self) -> PyResult<Vec<PyIon>> {
        let mut ions = Vec::new();
        let mut cm = self.cumulative_mass;

        for idx in 0..self.peptide.inner.sequence.len() - 1 {
            let r = self.peptide.inner.sequence[idx];
            let m = self.peptide.inner.modifications.get(idx).unwrap_or(&0.0);

            cm += match self.kind.inner {
                Kind::A | Kind::B | Kind::C => monoisotopic(r) + m,
                Kind::X | Kind::Y | Kind::Z => -(monoisotopic(r) + m),
            };

            ions.push(PyIon {
                inner: Ion {
                    kind: self.kind.inner.clone(),
                    monoisotopic_mass: cm,
                },
            });
        }
        Ok(ions)
    }
}

#[pymodule]
pub fn ion_series(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyKind>()?;
    m.add_class::<PyIon>()?;
    m.add_class::<PyIonSeries>()?;
    Ok(())
}
