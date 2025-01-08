use std::collections::HashMap;
use pyo3::prelude::*;
use unimod::unimod::{unimod_modifications_mass_numerical, unimod_modifications_mass, quantized_mass_to_unimod, quanzie_mass, title_to_unimod_id, modification_atomic_composition};

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

#[pyfunction]
fn title_to_unimod_ids() -> HashMap<&'static str, &'static str> {
    title_to_unimod_id()
}

#[pyfunction]
fn modification_atomic_compositions() -> HashMap<String, HashMap<&'static str, i32>> {
    modification_atomic_composition()
}

#[pymodule]
pub fn py_unimod(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(unimod_modification_to_mass_numerical, m)?)?;
    m.add_function(wrap_pyfunction!(unimod_modification_to_mass, m)?)?;
    m.add_function(wrap_pyfunction!(quantized_mass_to_unimod_candidates, m)?)?;
    m.add_function(wrap_pyfunction!(quanzied_mass, m)?)?;
    m.add_function(wrap_pyfunction!(title_to_unimod_ids, m)?)?;
    m.add_function(wrap_pyfunction!(modification_atomic_compositions, m)?)?;
    Ok(())
}