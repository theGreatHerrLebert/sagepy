use pyo3::prelude::*;

pub mod py_database;
pub mod py_enzyme;
pub mod py_fasta;
mod py_ion_series;
pub mod py_mass;
mod py_modification;
mod py_peptide;
mod py_scoring;
mod py_spectrum;
mod py_fdr;
mod py_lfq;
mod py_tmt;

use py_enzyme::enzyme;
use py_fasta::fasta;
use py_ion_series::ion_series;
use py_mass::mass;
use py_modification::modification;
use py_peptide::peptide;
use py_scoring::scoring;
use py_spectrum::spectrum;
use py_fdr::fdr;
use py_lfq::lfq;
use py_tmt::tmt;

#[pymodule]
fn sagepy_connector(py: Python, m: &PyModule) -> PyResult<()> {
    // py_mass submodule //
    let py_mass_submodule = PyModule::new(py, "py_mass")?;
    mass(py, &py_mass_submodule)?;
    m.add_submodule(py_mass_submodule)?;

    // py_enzyme submodule //
    let py_enzyme_submodule = PyModule::new(py, "py_enzyme")?;
    enzyme(py, &py_enzyme_submodule)?;
    m.add_submodule(py_enzyme_submodule)?;

    // py_fasta submodule //
    let py_fasta_submodule = PyModule::new(py, "py_fasta")?;
    fasta(py, &py_fasta_submodule)?;
    m.add_submodule(py_fasta_submodule)?;

    // py_peptide submodule //
    let py_peptide_submodule = PyModule::new(py, "py_peptide")?;
    peptide(py, &py_peptide_submodule)?;
    m.add_submodule(py_peptide_submodule)?;

    // py_ionseries submodule //
    let py_ionseries_submodule = PyModule::new(py, "py_ion_series")?;
    ion_series(py, &py_ionseries_submodule)?;
    m.add_submodule(py_ionseries_submodule)?;

    // py_modification submodule //
    let py_modification_submodule = PyModule::new(py, "py_modification")?;
    modification(py, &py_modification_submodule)?;
    m.add_submodule(py_modification_submodule)?;

    // py_database submodule //
    let py_database_submodule = PyModule::new(py, "py_database")?;
    py_database::database(py, &py_database_submodule)?;
    m.add_submodule(py_database_submodule)?;

    // py_spectrum submodule //
    let py_spectrum_submodule = PyModule::new(py, "py_spectrum")?;
    spectrum(py, &py_spectrum_submodule)?;
    m.add_submodule(py_spectrum_submodule)?;

    // py_scoring submodule //
    let py_scoring_submodule = PyModule::new(py, "py_scoring")?;
    scoring(py, &py_scoring_submodule)?;
    m.add_submodule(py_scoring_submodule)?;

    // py_fdr submodule //
    let py_fdr_submodule = PyModule::new(py, "py_fdr")?;
    fdr(py, &py_fdr_submodule)?;
    m.add_submodule(py_fdr_submodule)?;

    // py_lfq submodule //
    let py_lfq_submodule = PyModule::new(py, "py_lfq")?;
    lfq(py, &py_lfq_submodule)?;
    m.add_submodule(py_lfq_submodule)?;

    // py_tmt submodule //
    let py_tmt_submodule = PyModule::new(py, "py_tmt")?;
    tmt(py, &py_tmt_submodule)?;
    m.add_submodule(py_tmt_submodule)?;

    Ok(())
}
