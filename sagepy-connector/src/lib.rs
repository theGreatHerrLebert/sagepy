use pyo3::prelude::*;
use pyo3::wrap_pymodule;

pub mod py_database;
pub mod py_enzyme;
pub mod py_fasta;
pub mod py_ion_series;
pub mod py_mass;
pub mod py_modification;
pub mod py_peptide;
pub mod py_scoring;
pub mod py_spectrum;
pub mod py_fdr;
pub mod py_lfq;
pub mod py_tmt;
pub mod py_qfdr;
pub mod py_utility;
pub mod py_unimod;
pub mod utilities;
pub mod py_intensity;
pub mod py_retention_model;
pub mod py_retention_alignment;
pub mod py_mobility_model;
#[pymodule]
fn sagepy_connector(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {

    m.add_wrapped(wrap_pymodule!(py_mass::py_mass))?;
    m.add_wrapped(wrap_pymodule!(py_enzyme::py_enzyme))?;
    m.add_wrapped(wrap_pymodule!(py_fasta::py_fasta))?;
    m.add_wrapped(wrap_pymodule!(py_peptide::py_peptide))?;
    m.add_wrapped(wrap_pymodule!(py_ion_series::py_ion_series))?;
    m.add_wrapped(wrap_pymodule!(py_modification::py_modification))?;
    m.add_wrapped(wrap_pymodule!(py_database::py_database))?;
    m.add_wrapped(wrap_pymodule!(py_spectrum::py_spectrum))?;
    m.add_wrapped(wrap_pymodule!(py_scoring::py_scoring))?;
    m.add_wrapped(wrap_pymodule!(py_fdr::py_fdr))?;
    m.add_wrapped(wrap_pymodule!(py_lfq::py_lfq))?;
    m.add_wrapped(wrap_pymodule!(py_tmt::py_tmt))?;
    m.add_wrapped(wrap_pymodule!(py_qfdr::py_qfdr))?;
    m.add_wrapped(wrap_pymodule!(py_unimod::py_unimod))?;
    m.add_wrapped(wrap_pymodule!(py_utility::py_utility))?;
    m.add_wrapped(wrap_pymodule!(py_intensity::py_intensity))?;
    m.add_wrapped(wrap_pymodule!(py_retention_alignment::py_retention_alignment))?;

    Ok(())
}
