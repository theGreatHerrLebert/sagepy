use crate::py_peptide::PyPeptide;
use pyo3::prelude::*;
use sage_core::ion_series::{Ion, Kind};
use sage_core::mass::monoisotopic;

#[pyclass]
#[derive(Clone, Copy)]
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
        let seq_len = self.peptide.inner.sequence.len();
        let mut ions = Vec::with_capacity(seq_len.saturating_sub(1));
        let mut cm = self.cumulative_mass;

        for idx in 0..seq_len.saturating_sub(1) {
            let r = self.peptide.inner.sequence[idx];
            let m = self.peptide.inner.modifications.mass_at(idx);

            cm += match self.kind.inner {
                Kind::A | Kind::B | Kind::C => monoisotopic(r) + m,
                Kind::X | Kind::Y | Kind::Z => -(monoisotopic(r) + m),
            };

            ions.push(PyIon {
                inner: Ion {
                    kind: self.kind.inner,
                    monoisotopic_mass: cm,
                },
            });
        }
        Ok(ions)
    }
}

#[pymodule]
pub fn py_ion_series(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyKind>()?;
    m.add_class::<PyIon>()?;
    m.add_class::<PyIonSeries>()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use sage_core::enzyme::Digest;
    use sage_core::ion_series::IonSeries;
    use sage_core::peptide::Peptide;

    fn make_peptide(s: &str) -> Peptide {
        Peptide::try_from(Digest {
            sequence: s.into(),
            ..Default::default()
        })
        .unwrap()
    }

    #[test]
    fn test_kind_copy() {
        let k1 = PyKind { inner: Kind::B };
        let k2 = k1; // Copy
        assert_eq!(k1.kind_as_string(), k2.kind_as_string());
    }

    #[test]
    fn test_kind_from_string() {
        assert!(PyKind::new("b".to_string()).is_ok());
        assert!(PyKind::new("Y".to_string()).is_ok());
        assert!(PyKind::new("invalid".to_string()).is_err());
    }

    #[test]
    fn test_ion_series_length() {
        // Test directly with sage-core IonSeries
        let peptide = make_peptide("PEPTIDE");
        let ions: Vec<_> = IonSeries::new(&peptide, Kind::B).collect();

        // PEPTIDE has 7 residues, so 6 fragment ions (n-1)
        assert_eq!(ions.len(), 6);
    }

    #[test]
    fn test_b_ion_masses_increase() {
        let peptide = make_peptide("PEPTIDE");
        let ions: Vec<_> = IonSeries::new(&peptide, Kind::B).collect();

        // B ions should have increasing masses
        for i in 1..ions.len() {
            assert!(
                ions[i].monoisotopic_mass > ions[i - 1].monoisotopic_mass,
                "B ion masses should increase"
            );
        }
    }

    #[test]
    fn test_y_ion_masses_decrease() {
        let peptide = make_peptide("PEPTIDE");
        let ions: Vec<_> = IonSeries::new(&peptide, Kind::Y).collect();

        // Y ions should have decreasing masses (we iterate from N-term losing residues)
        for i in 1..ions.len() {
            assert!(
                ions[i].monoisotopic_mass < ions[i - 1].monoisotopic_mass,
                "Y ion masses should decrease"
            );
        }
    }
}
