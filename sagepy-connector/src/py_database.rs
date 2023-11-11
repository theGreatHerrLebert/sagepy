use numpy::{IntoPyArray, PyArray1};
use std::collections::HashMap;

use crate::py_enzyme::PyEnzymeParameters;
use crate::py_fasta::PyFasta;
use crate::py_ion_series::PyKind;
use crate::py_mass::PyTolerance;
use crate::py_modification::PyModificationSpecificity;
use crate::py_peptide::PyPeptide;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use sage_core::database::{
    Builder, EnzymeBuilder, IndexedDatabase, Parameters, PeptideIx, Theoretical,
};
use sage_core::fasta::Fasta;
use sage_core::ion_series::Kind;

#[pyclass]
#[derive(Clone)]
pub struct PyIndexedQuery {
    pub precursor_mass: f32,
    pub precursor_tolerance: PyTolerance,
    pub fragment_tolerance: PyTolerance,
    pub pre_idx_lo: usize,
    pub pre_idx_hi: usize,
}

#[pymethods]
impl PyIndexedQuery {
    #[new]
    pub fn new(
        precursor_mass: f32,
        precursor_tolerance: PyTolerance,
        fragment_tolerance: PyTolerance,
        pre_idx_lo: usize,
        pre_idx_hi: usize,
    ) -> PyResult<Self> {
        Ok(PyIndexedQuery {
            precursor_mass,
            precursor_tolerance,
            fragment_tolerance,
            pre_idx_lo,
            pre_idx_hi,
        })
    }

    #[getter]
    pub fn precursor_mass(&self) -> f32 {
        self.precursor_mass
    }

    #[getter]
    pub fn precursor_tolerance(&self) -> PyTolerance {
        self.precursor_tolerance.clone()
    }

    #[getter]
    pub fn fragment_tolerance(&self) -> PyTolerance {
        self.fragment_tolerance.clone()
    }

    #[getter]
    pub fn pre_idx_lo(&self) -> usize {
        self.pre_idx_lo
    }

    #[getter]
    pub fn pre_idx_hi(&self) -> usize {
        self.pre_idx_hi
    }
}

#[pyclass]
pub struct PyIndexedDatabase {
    pub inner: IndexedDatabase,
}

#[pymethods]
impl PyIndexedDatabase {
    #[new]
    pub fn new(
        peptides: Vec<PyPeptide>,
        fragments: Vec<PyTheoretical>,
        ion_kinds: Vec<PyKind>,
        min_value: Vec<f32>,
        potential_mods: Vec<(PyModificationSpecificity, f32)>,
        bucket_size: usize,
        generate_decoys: bool,
        decoy_tag: String,
    ) -> PyResult<Self> {
        Ok(PyIndexedDatabase {
            inner: IndexedDatabase {
                peptides: peptides.into_iter().map(|p| p.inner).collect(),
                fragments: fragments.into_iter().map(|f| f.inner).collect(),
                ion_kinds: ion_kinds.into_iter().map(|k| k.inner).collect(),
                min_value,
                potential_mods: potential_mods
                    .into_iter()
                    .map(|(k, v)| (k.inner, v))
                    .collect(),
                bucket_size,
                generate_decoys,
                decoy_tag,
            },
        })
    }

    #[staticmethod]
    pub fn from_parameters(parameters: PyParameters, fasta: PyFasta) -> PyResult<Self> {
        Ok(PyIndexedDatabase {
            inner: parameters.inner.build(fasta.inner),
        })
    }

    pub fn query(
        &self,
        precursor_mass: f32,
        precursor_tolerance: PyTolerance,
        fragment_tolerance: PyTolerance,
    ) -> PyResult<PyIndexedQuery> {
        let precursor_tolerance = precursor_tolerance.inner;
        let fragment_tolerance = fragment_tolerance.inner;
        let query = self
            .inner
            .query(precursor_mass, precursor_tolerance, fragment_tolerance);
        Ok(PyIndexedQuery {
            precursor_mass,
            precursor_tolerance: PyTolerance {
                inner: precursor_tolerance,
            },
            fragment_tolerance: PyTolerance {
                inner: fragment_tolerance,
            },
            pre_idx_lo: query.pre_idx_lo,
            pre_idx_hi: query.pre_idx_hi,
        })
    }

    #[getter]
    pub fn peptides(&self) -> Vec<PyPeptide> {
        self.inner
            .peptides
            .iter()
            .map(|p| PyPeptide { inner: p.clone() })
            .collect()
    }

    pub fn peptides_as_string(&self, _py: Python) -> Vec<String> {
        let peptides = self
            .inner
            .peptides
            .iter()
            .map(|p| String::from_utf8(p.sequence.clone().to_vec()))
            .collect::<Result<Vec<String>, _>>()
            .unwrap();
        peptides
    }

    pub fn mono_masses(&self, py: Python) -> Py<PyArray1<f32>> {
        let masses = self
            .inner
            .peptides
            .iter()
            .map(|p| p.monoisotopic)
            .collect::<Vec<f32>>();
        masses.into_pyarray(py).to_owned()
    }

    pub fn modifications(&self) -> Vec<(usize, Vec<f32>)> {
        let mut mods: Vec<(usize, Vec<f32>)> = Vec::new();

        for (i, p) in self.inner.peptides.iter().enumerate() {
            // only add if there are modifications based on non-zero entries in the modifications vector
            if p.modifications.iter().any(|m| *m != 0.0) {
                mods.push((i, p.modifications.clone()));
            }
        }
        mods
    }

    #[getter]
    pub fn num_peptides(&self) -> usize {
        self.inner.peptides.len()
    }

    #[getter]
    pub fn fragments(&self) -> Vec<PyTheoretical> {
        self.inner
            .fragments
            .iter()
            .map(|f| PyTheoretical { inner: f.clone() })
            .collect()
    }

    #[getter]
    pub fn num_fragments(&self) -> usize {
        self.inner.fragments.len()
    }

    #[getter]
    pub fn ion_kinds(&self) -> Vec<PyKind> {
        self.inner
            .ion_kinds
            .iter()
            .map(|k| PyKind { inner: k.clone() })
            .collect()
    }

    #[getter]
    pub fn min_value(&self) -> Vec<f32> {
        self.inner.min_value.clone()
    }

    #[getter]
    pub fn potential_mods(&self) -> Vec<(PyModificationSpecificity, f32)> {
        self.inner
            .potential_mods
            .iter()
            .map(|(k, v)| (PyModificationSpecificity { inner: k.clone() }, *v))
            .collect()
    }

    #[getter]
    pub fn bucket_size(&self) -> usize {
        self.inner.bucket_size
    }

    #[getter]
    pub fn generate_decoys(&self) -> bool {
        self.inner.generate_decoys
    }

    #[getter]
    pub fn decoy_tag(&self) -> String {
        self.inner.decoy_tag.clone()
    }
}

#[pyclass]
#[derive(Clone, Debug)]
pub struct PyEnzymeBuilder {
    pub inner: EnzymeBuilder,
}

#[pymethods]
impl PyEnzymeBuilder {
    #[new]
    pub fn new(
        missed_cleavages: Option<u8>,
        min_len: Option<usize>,
        max_len: Option<usize>,
        cleave_at: Option<String>,
        restrict: Option<char>,
        c_terminal: Option<bool>,
        semi_enzymatic: Option<bool>,
    ) -> PyResult<Self> {
        Ok(PyEnzymeBuilder {
            inner: EnzymeBuilder {
                missed_cleavages,
                min_len,
                max_len,
                cleave_at,
                restrict,
                c_terminal,
                semi_enzymatic,
            },
        })
    }

    #[staticmethod]
    pub fn from_default_trypsin() -> PyResult<Self> {
        Ok(PyEnzymeBuilder {
            inner: EnzymeBuilder::default(),
        })
    }

    pub fn get_enzyme_parameters(&self) -> PyResult<PyEnzymeParameters> {
        Ok(PyEnzymeParameters {
            inner: self.clone().inner.into(),
        })
    }

    #[getter]
    pub fn missed_cleavages(&self) -> Option<u8> {
        self.inner.missed_cleavages
    }

    #[getter]
    pub fn min_len(&self) -> Option<usize> {
        self.inner.min_len
    }

    #[getter]
    pub fn max_len(&self) -> Option<usize> {
        self.inner.max_len
    }

    #[getter]
    pub fn cleave_at(&self) -> Option<String> {
        self.inner.cleave_at.clone()
    }

    #[getter]
    pub fn restrict(&self) -> Option<char> {
        self.inner.restrict
    }

    #[getter]
    pub fn c_terminal(&self) -> Option<bool> {
        self.inner.c_terminal
    }
}

#[pyclass]
#[derive(Clone, Debug)]
pub struct PyPeptideIx {
    pub inner: PeptideIx,
}

#[pymethods]
impl PyPeptideIx {
    #[new]
    pub fn new(idx: u32) -> PyResult<Self> {
        Ok(PyPeptideIx {
            inner: PeptideIx(idx),
        })
    }

    #[getter]
    pub fn idx(&self) -> u32 {
        self.inner.0
    }
}

#[pyclass]
#[derive(Clone, Debug)]
pub struct PyTheoretical {
    pub inner: Theoretical,
}

#[pymethods]
impl PyTheoretical {
    #[new]
    pub fn new(idx: u32, fragment_mz: f32) -> PyResult<Self> {
        Ok(PyTheoretical {
            inner: Theoretical {
                peptide_index: PeptideIx(idx),
                fragment_mz,
            },
        })
    }

    #[getter]
    pub fn idx(&self) -> PyPeptideIx {
        PyPeptideIx {
            inner: self.inner.peptide_index,
        }
    }

    #[getter]
    pub fn fragment_mz(&self) -> f32 {
        self.inner.fragment_mz
    }
}

#[pyclass]
#[derive(Clone, Debug)]
pub struct PyParameters {
    pub inner: Parameters,
}

#[pymethods]
impl PyParameters {
    #[new]
    pub fn new(
        bucket_size: usize,
        py_enzyme_builder: PyEnzymeBuilder,
        fragment_min_mz: f32,
        fragment_max_mz: f32,
        peptide_min_mass: f32,
        peptide_max_mass: f32,
        min_ion_index: usize,
        static_mods: &PyDict,
        variable_mods: &PyDict,
        max_variable_mods: usize,
        decoy_tag: String,
        generate_decoys: bool,
        fasta: String,
        ion_kinds: Option<Vec<PyKind>>,
    ) -> PyResult<Self> {
        Ok(PyParameters {
            inner: Parameters {
                bucket_size,
                enzyme: py_enzyme_builder.inner,
                fragment_min_mz,
                fragment_max_mz,
                peptide_min_mass,
                peptide_max_mass,
                ion_kinds: ion_kinds
                    .map(|k| k.iter().map(|k| k.inner.clone()).collect())
                    .unwrap_or_else(|| vec![Kind::B, Kind::Y]),
                min_ion_index,
                static_mods: static_mods
                    .extract::<HashMap<PyModificationSpecificity, f32>>()
                    .unwrap()
                    .iter()
                    .map(|(k, v)| (k.inner.clone(), *v))
                    .collect(),
                variable_mods: variable_mods
                    .extract::<HashMap<PyModificationSpecificity, Vec<f32>>>()
                    .unwrap()
                    .iter()
                    .map(|(k, v)| (k.inner.clone(), v.clone()))
                    .collect(),
                max_variable_mods,
                decoy_tag,
                generate_decoys,
                fasta,
            },
        })
    }
    #[staticmethod]
    pub fn from_default() -> PyResult<Self> {
        Ok(PyParameters {
            inner: Builder::default().make_parameters(),
        })
    }

    pub fn digest(&self) -> PyResult<Vec<PyPeptide>> {
        let fasta = Fasta::parse(
            self.inner.fasta.clone(),
            self.inner.decoy_tag.clone(),
            self.inner.generate_decoys,
        );
        let digest = self.inner.digest(&fasta);
        Ok(digest.into_iter().map(|t| PyPeptide { inner: t }).collect())
    }

    pub fn build_indexed_database(&self) -> PyResult<PyIndexedDatabase> {
        Ok(PyIndexedDatabase {
            inner: self.inner.clone().build(Fasta::parse(
                self.inner.fasta.clone(),
                self.inner.decoy_tag.clone(),
                self.inner.generate_decoys,
            )),
        })
    }

    #[getter]
    pub fn bucket_size(&self) -> usize {
        self.inner.bucket_size
    }

    #[getter]
    pub fn enzyme_builder(&self) -> PyEnzymeBuilder {
        PyEnzymeBuilder {
            inner: self.inner.enzyme.clone(),
        }
    }

    #[getter]
    pub fn fragment_min_mz(&self) -> f32 {
        self.inner.fragment_min_mz
    }

    #[getter]
    pub fn fragment_max_mz(&self) -> f32 {
        self.inner.fragment_max_mz
    }

    #[getter]
    pub fn peptide_min_mass(&self) -> f32 {
        self.inner.peptide_min_mass
    }

    #[getter]
    pub fn peptide_max_mass(&self) -> f32 {
        self.inner.peptide_max_mass
    }

    #[getter]
    pub fn ion_kinds(&self) -> Vec<PyKind> {
        self.inner
            .ion_kinds
            .iter()
            .map(|k| PyKind { inner: k.clone() })
            .collect()
    }

    #[getter]
    pub fn min_ion_index(&self) -> usize {
        self.inner.min_ion_index
    }

    #[getter]
    pub fn static_mods(&self) -> HashMap<PyModificationSpecificity, f32> {
        self.inner
            .static_mods
            .iter()
            .map(|(k, v)| (PyModificationSpecificity { inner: k.clone() }, *v))
            .collect::<HashMap<PyModificationSpecificity, f32>>()
    }

    #[getter]
    pub fn variable_mods(&self) -> HashMap<PyModificationSpecificity, Vec<f32>> {
        self.inner
            .variable_mods
            .iter()
            .map(|(k, v)| (PyModificationSpecificity { inner: k.clone() }, v.clone()))
            .collect::<HashMap<PyModificationSpecificity, Vec<f32>>>()
    }

    #[getter]
    pub fn max_variable_mods(&self) -> usize {
        self.inner.max_variable_mods
    }

    #[getter]
    pub fn decoy_tag(&self) -> String {
        self.inner.decoy_tag.clone()
    }

    #[getter]
    pub fn generate_decoys(&self) -> bool {
        self.inner.generate_decoys
    }

    #[getter]
    pub fn fasta(&self) -> String {
        self.inner.fasta.clone()
    }
}

#[pymodule]
pub fn database(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyPeptideIx>()?;
    m.add_class::<PyTheoretical>()?;
    m.add_class::<PyParameters>()?;
    m.add_class::<PyEnzymeBuilder>()?;
    m.add_class::<PyIndexedDatabase>()?;
    m.add_class::<PyIndexedQuery>()?;
    Ok(())
}
