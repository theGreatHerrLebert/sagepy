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
use sage_core::peptide::Peptide;
use sage_core::scoring::{Feature, Scorer};
use crate::py_scoring::{make_sage_scorer, remove_duplicates, PyPsm, PyScorer};
use crate::py_spectrum::PyProcessedSpectrum;

use std::collections::BTreeMap;
use qfdrust::psm::Psm;
use rayon::prelude::*;
use rayon::ThreadPool;
use crate::utilities::sage_sequence_to_unimod_sequence;

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

    pub fn __getitem__(&self, index: PyPeptideIx) -> PyPeptide {
        PyPeptide {
            inner: self.inner[index.inner].clone(),
        }
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
        masses.into_pyarray(py).unbind()
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
    pub fn fragment_indices(&self, py: Python) -> Py<PyArray1<u32>> {
        let data: Vec<_> = self.inner
            .fragments
            .iter()
            .map(|f| f.peptide_index.0)
            .collect();
        data.into_pyarray(py).unbind()
    }

    #[getter]
    pub fn fragment_mzs(&self, py: Python) -> Py<PyArray1<f32>> {
        let data: Vec<_> = self.inner
            .fragments
            .iter()
            .map(|f| f.fragment_mz)
            .collect();
        data.into_pyarray(py).unbind()
    }

    pub fn fragment_dict(&self) -> HashMap<u32, Vec<f32>> {
        let mut fragment_dict: HashMap<u32, Vec<f32>> = HashMap::new();

        for fragment in &self.inner.fragments {
            fragment_dict.entry(fragment.peptide_index.0).or_insert(Vec::new()).push(fragment.fragment_mz);
        }

        // sort the fragment ions by m/z
        fragment_dict.iter_mut().for_each(|(_, v)| v.sort_by(|a, b| a.partial_cmp(b).unwrap()));

        fragment_dict
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
    #[pyo3(signature = (missed_cleavages=None, min_len=None, max_len=None, cleave_at=None, restrict=None, c_terminal=None, semi_enzymatic=None))]
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

    #[getter]
    pub fn semi_enzymatic(&self) -> Option<bool> {
        self.inner.semi_enzymatic
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
    #[pyo3(signature = (bucket_size, py_enzyme_builder, peptide_min_mass, peptide_max_mass, min_ion_index, static_mods, variable_mods, max_variable_mods, decoy_tag, generate_decoys, fasta, prefilter=None, prefilter_chunk_size=None, prefilter_low_memory=None, ion_kinds=None, _shuffle_decoys=None, _keep_ends=None))]
    pub fn new(
        bucket_size: usize,
        py_enzyme_builder: PyEnzymeBuilder,
        peptide_min_mass: f32,
        peptide_max_mass: f32,
        min_ion_index: usize,
        static_mods: &Bound<'_, PyDict>,
        variable_mods: &Bound<'_, PyDict>,
        max_variable_mods: usize,
        decoy_tag: String,
        generate_decoys: bool,
        fasta: String,
        prefilter: Option<bool>,
        prefilter_chunk_size: Option<usize>,
        prefilter_low_memory: Option<bool>,
        ion_kinds: Option<Vec<PyKind>>,
        _shuffle_decoys: Option<bool>,
        _keep_ends: Option<bool>,
    ) -> PyResult<Self> {
        Ok(PyParameters {
            inner: Parameters {
                bucket_size,
                enzyme: py_enzyme_builder.inner,
                peptide_min_mass,
                peptide_max_mass,
                ion_kinds: ion_kinds
                    .map(|k| k.iter().map(|k| k.inner.clone()).collect())
                    .unwrap_or_else(|| vec![Kind::B, Kind::Y]),
                min_ion_index,
                static_mods: static_mods
                    .extract::<HashMap<PyModificationSpecificity, f32>>()?
                    .iter()
                    .map(|(k, v)| (k.inner.clone(), *v))
                    .collect(),
                variable_mods: variable_mods
                    .extract::<HashMap<PyModificationSpecificity, Vec<f32>>>()?
                    .iter()
                    .map(|(k, v)| (k.inner.clone(), v.clone()))
                    .collect(),
                max_variable_mods,
                decoy_tag,
                generate_decoys,
                fasta,
                // shuffle_decoys: shuffle_decoys.unwrap_or(false),
                // keep_ends: keep_ends.unwrap_or(true),
                prefilter_chunk_size: prefilter_chunk_size.unwrap_or(0),
                prefilter: prefilter.unwrap_or(false),
                prefilter_low_memory: prefilter_low_memory.unwrap_or(false),
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
        let fasta = Fasta::parse(
            self.inner.fasta.clone(),
            self.inner.decoy_tag.clone(),
            self.inner.generate_decoys,
        );

        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            self.inner.clone().build(fasta)
        }));

        match result {
            Ok(db) => Ok(PyIndexedDatabase { inner: db }),
            Err(_) => Err(pyo3::exceptions::PyRuntimeError::new_err(
                "Rust panic occurred during indexed database generation. Check if input FASTA is empty or digestion parameters are too restrictive.",
            )),
        }
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

    pub fn prefilter_build_and_search_psm(
        &self,
        py: Python,
        spectra: Vec<PyProcessedSpectrum>,
        scorer_cfg: &PyScorer,
        chunk_size: usize,
        low_memory: bool,
        num_threads: usize,
        max_hits: usize,
    ) -> PyResult<(PyIndexedDatabase, BTreeMap<String, Vec<PyPsm>>, usize)> {
        py.allow_threads(|| -> PyResult<(PyIndexedDatabase, BTreeMap<String, Vec<PyPsm>>, usize)> {
            use rayon::prelude::*;
            use rayon::ThreadPoolBuilder;
            use std::sync::atomic::{AtomicBool, Ordering};

            // 0) Parse FASTA once
            let fasta = Fasta::parse(
                self.inner.fasta.clone(),
                self.inner.decoy_tag.clone(),
                self.inner.generate_decoys,
            );

            // 1) Thread pool
            let pool = ThreadPoolBuilder::new()
                .num_threads(num_threads)
                .build()
                .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

            // 2) Prefilter peptides
            let params = self.inner.clone();
            let mut kept_all: Vec<Peptide> = Vec::new();

            for fasta_chunk in fasta.iter_chunks(chunk_size) {
                let mut db = params.clone().build(fasta_chunk);

                // one AtomicBool per peptide in this chunk
                let keep: Vec<AtomicBool> = (0..db.peptides.len())
                    .map(|_| AtomicBool::new(false))
                    .collect();

                // scorer borrows chunk-db
                let sage_scorer = make_sage_scorer(scorer_cfg, &db);

                pool.install(|| {
                    spectra
                        .par_iter()
                        .filter(|s| s.inner.level == 2)
                        .for_each(|s| {
                            let _ = sage_scorer.quick_score(&s.inner, low_memory, keep.as_slice());
                        });
                });

                // drain kept peptides
                let kept = db
                    .peptides
                    .drain(..)
                    .enumerate()
                    .filter_map(|(ix, pep)| {
                        keep[ix].load(Ordering::Relaxed).then_some(pep)
                    })
                    .collect::<Vec<_>>();

                kept_all.extend(kept);
            }

            Parameters::reorder_peptides(&mut kept_all);

            let num_kept_peptides = kept_all.len();

            // 3) Build final DB once (moved out at end)
            let final_db = params.build_from_peptides(kept_all);

            // 4) Build PSM map using wrapper spectra (so CE is available)
            let psm_map: BTreeMap<String, Vec<PyPsm>> = build_psm_map_from_py_spectra(
                &final_db,
                scorer_cfg,
                spectra.as_slice(),
                &pool,
                max_hits,
            );

            Ok((PyIndexedDatabase { inner: final_db }, psm_map, num_kept_peptides))
        })
    }
}

pub(crate) fn build_psm_map_from_py_spectra(
    db: &IndexedDatabase,
    scorer_cfg: &PyScorer,
    spectra: &[PyProcessedSpectrum],
    pool: &ThreadPool,
    max_hits: usize,
) -> BTreeMap<String, Vec<PyPsm>> {
    // scorer borrows `db` (no clone)
    let scorer: Scorer<'_> = make_sage_scorer(scorer_cfg, db);

    // 1) Score all spectra (parallel)
    let features_per_spec: Vec<Vec<Feature>> = pool.install(|| {
        spectra
            .par_iter()
            .map(|s| scorer.score(&s.inner))
            .collect()
    });

    // 2) Build raw PSM map (parallel zip)
    let psm_map_raw: BTreeMap<String, Vec<Psm>> = pool.install(|| {
        spectra
            .par_iter()
            .zip(features_per_spec.into_par_iter())
            .map(|(spectrum, features)| {
                // MS1 intensity from inner (as in your code)
                let intensity_ms1: f32 = spectrum
                    .inner
                    .precursors
                    .iter()
                    .map(|p| p.intensity.unwrap_or(0.0))
                    .sum();

                // CE from OUTER wrapper (your correction)
                let collision_energy: f32 = spectrum
                    .collision_energies
                    .first()
                    .copied()
                    .flatten()
                    .unwrap_or(0.0f32);

                let mut psms: Vec<Psm> = Vec::with_capacity(features.len());

                for feature in &features {
                    let peptide = db[feature.peptide_idx].clone();

                    let intensity_ms2: f32 = feature.ms2_intensity;

                    let proteins: Vec<String> = peptide
                        .proteins
                        .iter()
                        .map(|arc| arc.as_ref().to_string())
                        .collect();

                    let sequence = std::str::from_utf8(&peptide.sequence).unwrap().to_string();
                    let sequence_with_mods = sage_sequence_to_unimod_sequence(
                        sequence.clone(),
                        &peptide.modifications,
                        &scorer_cfg.expected_mods,
                    );

                    // avoid calling reverse(true) twice
                    let decoy = peptide.reverse(true);
                    let sequence_decoy = std::str::from_utf8(&decoy.sequence).unwrap().to_string();
                    let sequence_decoy_with_mods = sage_sequence_to_unimod_sequence(
                        sequence_decoy.clone(),
                        &decoy.modifications,
                        &scorer_cfg.expected_mods,
                    );

                    let psm = Psm::new(
                        spectrum.inner.id.clone(),
                        feature.peptide_idx.0,
                        proteins,
                        feature.clone(),
                        Some(sequence),
                        Some(sequence_with_mods),
                        Some(sequence_decoy),
                        Some(sequence_decoy_with_mods),
                        Some(intensity_ms1),
                        Some(intensity_ms2),
                        Some(collision_energy),
                        None, // collision_energy_calibrated
                        None, // retention_time_projected
                        None, // prosit_predicted_intensities
                        None, // re_score
                    );

                    psms.push(psm);
                }

                (spectrum.inner.id.clone(), psms)
            })
            .collect()
    });

    // 3) Wrap into PyPsm
    let mut out: BTreeMap<String, Vec<PyPsm>> = BTreeMap::new();
    for (spec_id, psms) in psm_map_raw {
        out.insert(spec_id, psms.into_iter().map(|p| PyPsm { inner: p }).collect());
    }

    // 4) Sort + clip before de-dup (optional)
    for (_, psms) in out.iter_mut() {
        // descending hyperscore is what your remove_duplicates assumes anyway
        psms.sort_by(|a, b| {
            b.inner
                .sage_feature
                .hyperscore
                .partial_cmp(&a.inner.sage_feature.hyperscore)
                .unwrap()
                .then_with(|| a.inner.peptide_idx.partial_cmp(&b.inner.peptide_idx).unwrap())
                .then_with(|| a.inner.sage_feature.label.partial_cmp(&b.inner.sage_feature.label).unwrap())
        });
        psms.truncate(max_hits);
    }

    // 5) De-dup + final clip
    let mut out = remove_duplicates(out);
    for (_, psms) in out.iter_mut() {
        psms.truncate(max_hits);
    }

    out
}


#[pymodule]
pub fn py_database(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyPeptideIx>()?;
    m.add_class::<PyTheoretical>()?;
    m.add_class::<PyParameters>()?;
    m.add_class::<PyEnzymeBuilder>()?;
    m.add_class::<PyIndexedDatabase>()?;
    m.add_class::<PyIndexedQuery>()?;
    Ok(())
}
