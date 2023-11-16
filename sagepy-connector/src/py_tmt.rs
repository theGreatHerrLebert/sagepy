use pyo3::prelude::*;
use sage_core::tmt::Isobaric;

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
}

#[pymodule]
pub fn tmt(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyIsobaric>()?;
    Ok(())
}
