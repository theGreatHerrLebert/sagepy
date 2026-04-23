use crate::py_scoring::{PyFragments, PyPsm};
use crate::utilities::sage_sequence_to_unimod_sequence;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use qfdrust::psm::{compress_psms, decompress_psms, Psm};
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use sage_core::ion_series::Kind;
use sage_core::scoring::Fragments;
use std::collections::{BTreeMap, HashMap, HashSet};

/// Converts a cosine similarity to an angle similarity.
/// The angle similarity is calculated as 1 - angle / pi.
///
/// # Arguments
///
/// * `cosim` - A f32 representing the cosine similarity.
///
/// # Returns
///
/// * A f32 representing the angle similarity.
///
#[pyfunction]
pub fn cosim_to_spectral_angle(cosim: f32) -> f32 {
    let angle = (1.0 - cosim).acos();
    1.0 - angle / std::f32::consts::PI
}

fn validate_equal_length<T, U>(left: &[T], right: &[U], left_name: &str, right_name: &str) -> PyResult<()> {
    if left.len() != right.len() {
        return Err(PyValueError::new_err(format!(
            "{left_name} and {right_name} must have the same length"
        )));
    }
    Ok(())
}

fn ppm_errors(measured_values: &[f64], reference_values: &[f64]) -> Vec<f64> {
    measured_values
        .iter()
        .zip(reference_values.iter())
        .map(|(measured, reference)| ((*measured - *reference) / *reference) * 1_000_000.0)
        .collect()
}

fn median(values: &mut [f64]) -> f64 {
    values.sort_unstable_by(f64::total_cmp);
    let mid = values.len() / 2;
    if values.len() % 2 == 0 {
        (values[mid - 1] + values[mid]) / 2.0
    } else {
        values[mid]
    }
}

fn std(sample: &[f64]) -> f64 {
    let mean = sample.iter().sum::<f64>() / sample.len() as f64;
    let variance = sample
        .iter()
        .map(|value| {
            let delta = *value - mean;
            delta * delta
        })
        .sum::<f64>()
        / sample.len() as f64;
    variance.sqrt()
}

fn kde_pdf(sample: &[f64], bandwidth: f64, x: f64) -> f64 {
    if sample.is_empty() || bandwidth <= 0.0 {
        return 0.0;
    }

    let constant = (2.0 * std::f64::consts::PI).sqrt() * bandwidth * sample.len() as f64;
    let sum_pdf = sample
        .par_iter()
        .map(|value| (-0.5 * ((x - *value) / bandwidth).powi(2)).exp())
        .sum::<f64>();
    sum_pdf / constant
}

fn posterior_error_inner(pep_bins: &[f64], min_score: f64, score_step: f64, score: f64) -> f64 {
    let last_index = pep_bins.len().saturating_sub(1) as isize;
    let mut bin_lo = ((score - min_score) / score_step) as isize;
    if bin_lo < 0 {
        bin_lo = 0;
    }
    if bin_lo > last_index {
        bin_lo = last_index;
    }

    let bin_lo = bin_lo as usize;
    let bin_hi = (bin_lo + 1).min(pep_bins.len() - 1);

    let lower = pep_bins[bin_lo];
    let upper = pep_bins[bin_hi];
    let bin_lo_score = bin_lo as f64 * score_step + min_score;
    let linear = (score - bin_lo_score) / score_step;
    lower + (upper - lower) * linear
}

#[pyfunction]
pub fn calculate_ppm_error(measured_value: f64, reference_value: f64) -> f64 {
    ((measured_value - reference_value) / reference_value) * 1_000_000.0
}

#[pyfunction]
pub fn calculate_ppms(measured_values: Vec<f64>, reference_values: Vec<f64>) -> PyResult<Vec<f64>> {
    validate_equal_length(
        measured_values.as_slice(),
        reference_values.as_slice(),
        "measured_values",
        "reference_values",
    )?;
    Ok(ppm_errors(measured_values.as_slice(), reference_values.as_slice()))
}

#[pyfunction]
pub fn mean_ppm(mz_observed: Vec<f64>, mz_calculated: Vec<f64>) -> PyResult<f64> {
    validate_equal_length(
        mz_observed.as_slice(),
        mz_calculated.as_slice(),
        "mz_observed",
        "mz_calculated",
    )?;
    if mz_observed.is_empty() {
        return Err(PyValueError::new_err("mz_observed cannot be empty"));
    }

    let errors = ppm_errors(mz_observed.as_slice(), mz_calculated.as_slice());
    Ok(errors.iter().sum::<f64>() / errors.len() as f64)
}

#[pyfunction]
pub fn median_ppm(mz_observed: Vec<f64>, mz_calculated: Vec<f64>) -> PyResult<f64> {
    validate_equal_length(
        mz_observed.as_slice(),
        mz_calculated.as_slice(),
        "mz_observed",
        "mz_calculated",
    )?;
    if mz_observed.is_empty() {
        return Err(PyValueError::new_err("mz_observed cannot be empty"));
    }

    let mut errors = ppm_errors(mz_observed.as_slice(), mz_calculated.as_slice());
    Ok(median(errors.as_mut_slice()))
}

#[pyfunction]
#[pyo3(signature = (scores, decoys, bins=1000, bw_adjust=1.0, monotonic=true))]
pub fn calculate_pep_single(
    scores: Vec<f64>,
    decoys: Vec<bool>,
    bins: usize,
    bw_adjust: f64,
    monotonic: bool,
) -> PyResult<(Vec<f64>, f64, f64)> {
    validate_equal_length(scores.as_slice(), decoys.as_slice(), "scores", "decoys")?;
    if scores.is_empty() {
        return Err(PyValueError::new_err("scores cannot be empty"));
    }
    if bins < 2 {
        return Err(PyValueError::new_err("bins must be at least 2"));
    }

    let decoy_scores: Vec<f64> = scores
        .iter()
        .zip(decoys.iter())
        .filter_map(|(score, is_decoy)| is_decoy.then_some(*score))
        .collect();
    let target_scores: Vec<f64> = scores
        .iter()
        .zip(decoys.iter())
        .filter_map(|(score, is_decoy)| (!*is_decoy).then_some(*score))
        .collect();

    if decoy_scores.is_empty() || target_scores.is_empty() {
        return Err(PyValueError::new_err(
            "scores must contain at least one decoy and one target",
        ));
    }

    let pi = decoy_scores.len() as f64 / scores.len() as f64;
    let sigma_d = std(decoy_scores.as_slice());
    let sigma_t = std(target_scores.as_slice());
    let bandwidth_d = bw_adjust * sigma_d * (4.0 / 3.0 / decoy_scores.len() as f64).powf(0.2);
    let bandwidth_t = bw_adjust * sigma_t * (4.0 / 3.0 / target_scores.len() as f64).powf(0.2);

    let min_score = scores.iter().copied().fold(f64::INFINITY, f64::min);
    let max_score = scores.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let score_step = (max_score - min_score) / (bins - 1) as f64;

    if score_step == 0.0 {
        return Err(PyValueError::new_err(
            "scores must span more than one distinct value",
        ));
    }

    let mut pep_bins: Vec<f64> = (0..bins)
        .into_par_iter()
        .map(|bin_idx| {
            let score = min_score + bin_idx as f64 * score_step;
            let decoy_pdf = kde_pdf(decoy_scores.as_slice(), bandwidth_d, score) * pi;
            let target_pdf = kde_pdf(target_scores.as_slice(), bandwidth_t, score) * (1.0 - pi);
            let denominator = decoy_pdf + target_pdf;
            if denominator == 0.0 {
                0.0
            } else {
                decoy_pdf / denominator
            }
        })
        .collect();

    if monotonic {
        for idx in (0..pep_bins.len() - 1).rev() {
            pep_bins[idx] = pep_bins[idx].max(pep_bins[idx + 1]);
        }
    }

    Ok((pep_bins, min_score, score_step))
}

#[pyfunction]
pub fn posterior_error(pep_bins: Vec<f64>, min_score: f64, score_step: f64, score: f64) -> PyResult<f64> {
    if pep_bins.is_empty() {
        return Err(PyValueError::new_err("pep_bins cannot be empty"));
    }
    if score_step == 0.0 {
        return Err(PyValueError::new_err("score_step cannot be zero"));
    }

    Ok(posterior_error_inner(
        pep_bins.as_slice(),
        min_score,
        score_step,
        score,
    ))
}

#[pyfunction]
#[pyo3(signature = (scores, decoys, bins=1000, bw_adjust=1.0, monotonic=true))]
pub fn calculate_pep(
    scores: Vec<f64>,
    decoys: Vec<bool>,
    bins: usize,
    bw_adjust: f64,
    monotonic: bool,
) -> PyResult<Vec<f64>> {
    let (pep_bins, min_score, score_step) = calculate_pep_single(
        scores.clone(),
        decoys,
        bins,
        bw_adjust,
        monotonic,
    )?;

    Ok(scores
        .into_par_iter()
        .map(|score| posterior_error_inner(pep_bins.as_slice(), min_score, score_step, score))
        .collect())
}

/// Reshape the flat prosit array into a 3D array of shape (29, 2, 3)
///
/// # Arguments
///
/// * `flat_array` - a vector of f64 representing the flat prosit array
///
/// # Returns
///
/// * `Vec<Vec<Vec<f64>>>` - a 3D array of shape (29, 2, 3)
///
pub fn reshape_prosit_array(flat_array: Vec<f32>) -> Vec<Vec<Vec<f32>>> {
    let mut array_return: Vec<Vec<Vec<f32>>> = vec![vec![vec![0.0; 3]; 2]; 29];
    let mut ptr = 0;

    for c in 0..3 {
        for row in 0..29 {
            // Fill in the Y ion values
            array_return[row][0][c] = flat_array[ptr];
            ptr += 1;
        }
        for row in 0..29 {
            // Fill in the B ion values
            array_return[row][1][c] = flat_array[ptr];
            ptr += 1;
        }
    }

    array_return
}

#[pyfunction]
pub fn flat_prosit_array_to_fragments_map(
    flat_intensities: Vec<f32>,
) -> BTreeMap<(u32, i32, i32), f32> {
    // Reshape the flat prosit array into a 3D array of shape (29, 2, 3)
    let reshaped_intensities = reshape_prosit_array(flat_intensities);

    // create hashmap of (kind, charge, ordinal) -> intensity
    let mut fragments: BTreeMap<(u32, i32, i32), f32> = BTreeMap::new();
    for z in 1..=3 {
        let intensity_b: Vec<f32> = reshaped_intensities[..]
            .iter()
            .map(|x| x[1][z as usize - 1])
            .collect();
        for i in 1..=29 {
            let intensity = intensity_b[i as usize - 1];
            if intensity >= 0.0 {
                fragments.insert((0, z, i), intensity);
            }
        }

        let intensity_y: Vec<f32> = reshaped_intensities[..]
            .iter()
            .map(|x| x[0][z as usize - 1])
            .collect();
        for i in 1..=29 {
            let intensity = intensity_y[i as usize - 1];
            if intensity >= 0.0 {
                fragments.insert((1, z, i), intensity);
            }
        }
    }
    fragments
}

#[pyfunction]
pub fn py_fragments_to_fragments_map(
    fragments: &PyFragments,
    normalize: bool,
) -> BTreeMap<(u32, i32, i32), f32> {
    let mut fragments_map: BTreeMap<(u32, i32, i32), f32> = BTreeMap::new();

    let max_intensity = fragments
        .inner
        .intensities
        .iter()
        .cloned()
        .fold(f32::NEG_INFINITY, f32::max);

    for i in 0..fragments.inner.mz_calculated.len() {
        let kind = match fragments.inner.kinds[i] {
            Kind::B => 0,
            Kind::Y => 1,
            _ => panic!("Invalid ion kind"),
        };

        let intensity = if normalize {
            fragments.inner.intensities[i] / max_intensity
        } else {
            fragments.inner.intensities[i]
        };

        fragments_map.insert(
            (
                kind,
                fragments.inner.charges[i],
                fragments.inner.fragment_ordinals[i],
            ),
            intensity,
        );
    }
    fragments_map
}

pub fn _map_to_py_fragments(
    fragments: &HashMap<(u32, i32, i32), f32>,
    mz_calculated: Vec<f32>,
    mz_experimental: Vec<f32>,
) -> PyFragments {
    let cap = fragments.len();
    let mut kinds: Vec<Kind> = Vec::with_capacity(cap);
    let mut ordinals: Vec<i32> = Vec::with_capacity(cap);
    let mut charges: Vec<i32> = Vec::with_capacity(cap);
    let mut intensities: Vec<f32> = Vec::with_capacity(cap);

    for ((kind_id, ordinal, charge), intensity) in fragments {
        let kind = match kind_id {
            0 => Kind::B,
            1 => Kind::Y,
            _ => panic!("Invalid ion kind"),
        };
        kinds.push(kind);
        ordinals.push(*ordinal);
        charges.push(*charge);
        intensities.push(*intensity);
    }

    PyFragments {
        inner: Fragments {
            mz_calculated,
            mz_experimental,
            kinds,
            fragment_ordinals: ordinals,
            charges,
            intensities,
        },
    }
}

#[pyfunction]
pub fn psms_to_json(psms: Vec<PyPsm>, num_threads: usize) -> Vec<String> {
    let thread_pool = ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();

    thread_pool.install(|| {
        psms.par_iter()
            .map(|psm| serde_json::to_string(&psm.inner).unwrap())
            .collect()
    })
}

#[pyfunction]
pub fn psms_to_json_bin(psms: Vec<PyPsm>) -> Vec<u8> {
    let inner_psms: Vec<Psm> = psms.into_iter().map(|psm| psm.inner).collect();
    bincode::serialize(&inner_psms).unwrap()
}

#[pyfunction]
pub fn json_bin_to_psms(json_bin: Vec<u8>) -> Vec<PyPsm> {
    let inner_psms: Vec<Psm> = bincode::deserialize(&json_bin).unwrap();
    inner_psms
        .into_iter()
        .map(|psm| PyPsm { inner: psm })
        .collect()
}

#[pyfunction]
pub fn sage_sequence_to_unimod(
    sequence: String,
    modifications: Vec<f32>,
    expected_modifications: HashSet<String>,
) -> String {
    sage_sequence_to_unimod_sequence(sequence, &modifications, &expected_modifications)
}

#[pyfunction]
pub fn psms_to_feature_matrix(psms: Vec<PyPsm>, num_threads: usize) -> Vec<Vec<f64>> {
    let thread_pool = ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();

    thread_pool.install(|| {
        psms.par_iter()
            .map(|psm| psm.inner.get_feature_vector())
            .collect()
    })
}

#[pyfunction]
pub fn get_psm_sequences_par(psms: Vec<PyPsm>, num_threads: usize) -> Vec<String> {
    let thread_pool = ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();

    thread_pool.install(|| {
        psms.par_iter()
            .map(|psm| psm.inner.sequence.as_ref().unwrap().sequence.clone())
            .collect()
    })
}

#[pyfunction]
pub fn get_psm_peptide_idx_par(psms: Vec<PyPsm>, num_threads: usize) -> Vec<u32> {
    let thread_pool = ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();

    thread_pool.install(|| {
        psms.par_iter()
            .map(|psm| psm.inner.sage_feature.peptide_idx.0)
            .collect()
    })
}

#[pyfunction]
pub fn get_psm_sequences_modified_par(psms: Vec<PyPsm>, num_threads: usize) -> Vec<String> {
    let thread_pool = ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();

    thread_pool.install(|| {
        psms.par_iter()
            .map(|psm| {
                psm.inner
                    .sequence_modified
                    .as_ref()
                    .unwrap()
                    .sequence
                    .clone()
            })
            .collect()
    })
}

#[pyfunction]
pub fn get_psm_sequences_decoy_par(psms: Vec<PyPsm>, num_threads: usize) -> Vec<String> {
    let thread_pool = ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();

    thread_pool.install(|| {
        psms.par_iter()
            .map(|psm| psm.inner.sequence_decoy.as_ref().unwrap().sequence.clone())
            .collect()
    })
}

#[pyfunction]
pub fn get_psm_sequences_decoy_modified_par(psms: Vec<PyPsm>, num_threads: usize) -> Vec<String> {
    let thread_pool = ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();

    thread_pool.install(|| {
        psms.par_iter()
            .map(|psm| {
                let sequence = match &psm.inner.sequence_decoy_modified {
                    Some(seq) => seq.sequence.clone(),
                    None => "".to_string(),
                };

                sequence
            })
            .collect()
    })
}

#[pyfunction]
pub fn get_psm_spec_idx_par(psms: Vec<PyPsm>, num_threads: usize) -> Vec<String> {
    let thread_pool = ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();

    thread_pool.install(|| {
        psms.par_iter()
            .map(|psm| psm.inner.spec_idx.clone())
            .collect()
    })
}

#[pyfunction]
pub fn get_psm_proteins_par(psms: Vec<PyPsm>, num_threads: usize) -> Vec<Vec<String>> {
    let thread_pool = ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();

    thread_pool.install(|| {
        psms.par_iter()
            .map(|psm| psm.inner.proteins.clone())
            .collect()
    })
}

#[pyfunction]
pub fn py_compress_psms(psms: Vec<PyPsm>) -> Vec<u8> {
    let inner_psms: Vec<Psm> = psms.into_iter().map(|psm| psm.inner).collect();
    compress_psms(&inner_psms).unwrap()
}

#[pyfunction]
pub fn py_decompress_psms(psms_bin: Vec<u8>) -> Vec<PyPsm> {
    let inner_psms: Vec<Psm> = decompress_psms(psms_bin.as_slice()).unwrap();
    inner_psms
        .into_iter()
        .map(|psm| PyPsm { inner: psm })
        .collect()
}

#[pymodule]
pub fn py_utility(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(calculate_ppm_error, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_ppms, m)?)?;
    m.add_function(wrap_pyfunction!(mean_ppm, m)?)?;
    m.add_function(wrap_pyfunction!(median_ppm, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_pep_single, m)?)?;
    m.add_function(wrap_pyfunction!(posterior_error, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_pep, m)?)?;
    m.add_function(wrap_pyfunction!(flat_prosit_array_to_fragments_map, m)?)?;
    m.add_function(wrap_pyfunction!(py_fragments_to_fragments_map, m)?)?;
    m.add_function(wrap_pyfunction!(psms_to_json, m)?)?;
    m.add_function(wrap_pyfunction!(psms_to_json_bin, m)?)?;
    m.add_function(wrap_pyfunction!(json_bin_to_psms, m)?)?;
    m.add_function(wrap_pyfunction!(cosim_to_spectral_angle, m)?)?;
    m.add_function(wrap_pyfunction!(sage_sequence_to_unimod, m)?)?;
    m.add_function(wrap_pyfunction!(psms_to_feature_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(get_psm_sequences_par, m)?)?;
    m.add_function(wrap_pyfunction!(get_psm_sequences_modified_par, m)?)?;
    m.add_function(wrap_pyfunction!(get_psm_sequences_decoy_par, m)?)?;
    m.add_function(wrap_pyfunction!(get_psm_sequences_decoy_modified_par, m)?)?;
    m.add_function(wrap_pyfunction!(get_psm_spec_idx_par, m)?)?;
    m.add_function(wrap_pyfunction!(get_psm_peptide_idx_par, m)?)?;
    m.add_function(wrap_pyfunction!(get_psm_proteins_par, m)?)?;
    m.add_function(wrap_pyfunction!(py_compress_psms, m)?)?;
    m.add_function(wrap_pyfunction!(py_decompress_psms, m)?)?;
    Ok(())
}
