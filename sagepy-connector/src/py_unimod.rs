use std::collections::HashMap;
use pyo3::prelude::*;
use unimod::unimod::{unimod_modifications_mass_numerical, unimod_modifications_mass, quantized_mass_to_unimod, quanzie_mass};

#[pyfunction]
fn unimod_modification_to_mass_numerical() -> HashMap<u32, f64> {
    unimod_modifications_mass_numerical()
}

#[pyfunction]
fn unimod_modification_to_mass() -> HashMap<&'static str, f64> {
    unimod_modifications_mass()
}

#[pyfunction]
fn quantized_mass_to_unimod_candidates() -> HashMap<i32, Vec<&'static str>> {
    quantized_mass_to_unimod()
}

#[pyfunction]
fn quanzied_mass(mass: f32) -> i32 {
    quanzie_mass(mass)
}

#[pymodule]
pub fn unimodifications(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(unimod_modification_to_mass_numerical, m)?)?;
    m.add_function(wrap_pyfunction!(unimod_modification_to_mass, m)?)?;
    m.add_function(wrap_pyfunction!(quantized_mass_to_unimod_candidates, m)?)?;
    m.add_function(wrap_pyfunction!(quanzied_mass, m)?)?;
    Ok(())
}