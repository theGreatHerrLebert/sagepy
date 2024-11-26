use std::collections::{BTreeMap, HashMap};
use std::f64::consts::LN_2;
use ndarray::Array1;
use ndarray::Zip;
use sage_core::ion_series::Kind;
use sage_core::scoring::Fragments;
use serde::{Deserialize, Serialize};

fn cosine_similarity(vec1: &Vec<f32>, vec2: &Vec<f32>, epsilon: f32) -> Option<f32> {
    if vec1.len() != vec2.len() || vec1.is_empty() {
        return None;
    }

    // filter the intensities based on the epsilon value
    let valid_ion_mask: Vec<bool> = vec2.iter().map(|&x| x > epsilon).collect();
    let vec1: Vec<f32> = vec1.iter().zip(&valid_ion_mask).filter_map(|(&x, &valid)| if valid { Some(x) } else { None }).collect();
    let vec2: Vec<f32> = vec2.iter().zip(&valid_ion_mask).filter_map(|(&x, &valid)| if valid { Some(x) } else { None }).collect();

    let dot_product: f32 = vec1.iter().zip(vec2.iter()).map(|(a, b)| a * b).sum();
    let magnitude_vec1: f32 = vec1.iter().map(|x| x.powi(2)).sum::<f32>().sqrt();
    let magnitude_vec2: f32 = vec2.iter().map(|x| x.powi(2)).sum::<f32>().sqrt();

    if magnitude_vec1 == 0.0 || magnitude_vec2 == 0.0 {
        return Some(0.0);
    }

    Some(dot_product / (magnitude_vec1 * magnitude_vec2))
}

fn cosim_to_spectral_angle(cosim: f32) -> f32 {
    let angle = cosim.acos(); // Use cosim directly
    1.0 - angle / std::f32::consts::PI
}


fn rank_ties(values: &Array1<f32>) -> Array1<f32> {
    let mut ranks = Array1::zeros(values.len());
    let mut value_to_indices: HashMap<i32, Vec<usize>> = HashMap::new();
    let quantization_factor = 1_000_000; // Adjust this factor as needed

    for (i, &value) in values.iter().enumerate() {
        let quantized_value = (value * quantization_factor as f32) as i32;
        value_to_indices.entry(quantized_value).or_default().push(i);
    }

    let mut sorted_values: Vec<_> = value_to_indices.keys().cloned().collect();
    sorted_values.sort();

    let mut rank = 1.0;
    for quantized_value in sorted_values {
        let indices = &value_to_indices[&quantized_value];
        let rank_sum: f32 = indices.iter().map(|_| rank).sum();
        let average_rank = rank_sum / indices.len() as f32;
        for &index in indices {
            ranks[index] = average_rank;
        }
        rank += indices.len() as f32;
    }

    ranks
}

fn pearson_correlation(observed_intensities: &[f32], predicted_intensities: &[f32], epsilon: f32) -> f32 {
    // Filter the intensities based on the epsilon value
    let valid_ion_mask: Vec<bool> = predicted_intensities.iter().map(|&x| x > epsilon).collect();
    let observed_filtered: Vec<f32> = observed_intensities.iter().zip(&valid_ion_mask).filter_map(|(&x, &valid)| if valid { Some(x) } else { None }).collect();
    let predicted_filtered: Vec<f32> = predicted_intensities.iter().zip(&valid_ion_mask).filter_map(|(&x, &valid)| if valid { Some(x) } else { None }).collect();

    // Remove NaN values
    let observed_filtered: Vec<f32> = observed_filtered.into_iter().filter(|&x| !x.is_nan()).collect();
    let predicted_filtered: Vec<f32> = predicted_filtered.into_iter().filter(|&x| !x.is_nan()).collect();

    if observed_filtered.len() <= 2 || predicted_filtered.len() <= 2 {
        return 0.0;
    }

    // Convert to ndarray
    let observed_array = Array1::from(observed_filtered);
    let predicted_array = Array1::from(predicted_filtered);

    // Calculate means
    let mean_observed = observed_array.mean().unwrap();
    let mean_predicted = predicted_array.mean().unwrap();

    // Calculate covariance and standard deviations
    let covariance = Zip::from(&observed_array)
        .and(&predicted_array)
        .map_collect(|&o, &p| (o - mean_observed) * (p - mean_predicted))
        .sum();
    let std_observed = (observed_array.mapv(|o| (o - mean_observed).powi(2)).sum()).sqrt();
    let std_predicted = (predicted_array.mapv(|p| (p - mean_predicted).powi(2)).sum()).sqrt();

    // Calculate Pearson correlation
    let corr = covariance / (std_observed * std_predicted);

    if corr.is_nan() {
        0.0
    } else {
        corr as f32
    }
}

fn spearman_correlation(observed_intensities: &[f32], predicted_intensities: &[f32], epsilon: f32) -> f32 {
    // Filter the intensities based on the epsilon value
    let valid_ion_mask: Vec<bool> = predicted_intensities.iter().map(|&x| x > epsilon).collect();
    let observed_filtered: Vec<f32> = observed_intensities.iter().zip(&valid_ion_mask).filter_map(|(&x, &valid)| if valid { Some(x) } else { None }).collect();
    let predicted_filtered: Vec<f32> = predicted_intensities.iter().zip(&valid_ion_mask).filter_map(|(&x, &valid)| if valid { Some(x) } else { None }).collect();

    // Remove NaN values
    let observed_filtered: Vec<f32> = observed_filtered.into_iter().filter(|&x| !x.is_nan()).collect();
    let predicted_filtered: Vec<f32> = predicted_filtered.into_iter().filter(|&x| !x.is_nan()).collect();

    if observed_filtered.len() <= 2 || predicted_filtered.len() <= 2 {
        return 0.0;
    }

    // Convert to ndarray
    let observed_array = Array1::from(observed_filtered);
    let predicted_array = Array1::from(predicted_filtered);

    // Rank the values
    let observed_ranks = rank_ties(&observed_array);
    let predicted_ranks = rank_ties(&predicted_array);

    // Calculate Pearson correlation of the ranks
    pearson_correlation(&observed_ranks.to_vec(), &predicted_ranks.to_vec(), epsilon)
}


fn spectral_entropy_similarity(observed_intensities: &Vec<f32>, predicted_intensities: &Vec<f32>, epsilon: f32) -> f32 {
    // Filter the intensities based on the epsilon value
    let valid_ion_mask: Vec<bool> = predicted_intensities.iter().map(|&x| x > epsilon).collect();
    let observed_filtered: Vec<f32> = observed_intensities.iter().zip(&valid_ion_mask).filter_map(|(&x, &valid)| if valid { Some(x) } else { None }).collect();
    let predicted_filtered: Vec<f32> = predicted_intensities.iter().zip(&valid_ion_mask).filter_map(|(&x, &valid)| if valid { Some(x) } else { None }).collect();

    // Calculate the entropy for the observed, predicted, and merged intensities
    let entropy_merged = calculate_entropy(&observed_filtered.iter().zip(&predicted_filtered).map(|(&o, &p)| o + p).collect::<Vec<f32>>());
    let entropy_obs = calculate_entropy(&observed_filtered);
    let entropy_pred = calculate_entropy(&predicted_filtered);

    // Calculate the spectral entropy similarity
    let entropy = 1.0 - (2.0 * entropy_merged - entropy_obs - entropy_pred) / (2.0 * LN_2 as f32);

    // Handle cases where the entropy is NaN (set them to 0)
    if entropy.is_nan() {
        0.0
    } else {
        entropy
    }
}

fn calculate_entropy(intensities: &Vec<f32>) -> f32 {
    let sum: f32 = intensities.iter().sum();
    if sum == 0.0 {
        return 0.0;
    }

    intensities.iter().map(|&x| {
        let p = x / sum;
        if p > 0.0 {
            -p * p.ln()
        } else {
            0.0
        }
    }).sum()
}

pub fn flat_prosit_array_to_fragments_map(flat_intensities: Vec<f32>) -> BTreeMap<(u32, i32, i32), f32> {
    // Reshape the flat prosit array into a 3D array of shape (29, 2, 3)
    let reshaped_intensities = reshape_prosit_array(flat_intensities);

    // create hashmap of (kind, charge, ordinal) -> intensity
    let mut fragments: BTreeMap<(u32, i32, i32), f32> = BTreeMap::new();
    for z in 1..=3 {
        let intensity_b: Vec<f32> = reshaped_intensities[..].iter().map(|x| x[1][z as usize - 1]).collect();
        for i in 1..=29 {
            let intensity = intensity_b[i as usize - 1];
            if intensity >= 0.0 {
                fragments.insert((0, z, i), intensity);
            }
        }

        let intensity_y: Vec<f32> = reshaped_intensities[..].iter().map(|x| x[0][z as usize - 1]).collect();
        for i in 1..=29 {
            let intensity = intensity_y[i as usize - 1];
            if intensity >= 0.0 {
                fragments.insert((1, z, i), intensity);
            }
        }
    }
    fragments
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

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FragmentIntensityPrediction {
    pub fragments: Fragments,
    pub prosit_intensity_predicted: Vec<f32>,
}

impl FragmentIntensityPrediction {
    pub fn new(
        fragments: Fragments,
        prosit_intensity_predicted:Vec<f32>,
    ) -> Self {
        FragmentIntensityPrediction {
            fragments,
            prosit_intensity_predicted,
        }
    }
    pub fn prosit_intensity_to_fragments_map(&self) -> BTreeMap<(u32, i32, i32), f32> {
        flat_prosit_array_to_fragments_map(self.prosit_intensity_predicted.clone())
    }

    pub fn observed_intensity_to_fragments_map(&self) -> BTreeMap<(u32, i32, i32), f32> {
        let mut fragments: BTreeMap<(u32, i32, i32), f32> = BTreeMap::new();
        let max_intensity = self.fragments.intensities.iter().cloned().fold(f32::NEG_INFINITY, f32::max);

        for i in 0..self.fragments.mz_calculated.len() {
            let kind = match self.fragments.kinds[i] {
                Kind::B => 0,
                Kind::Y => 1,
                _ => panic!("Invalid kind"),
            };
            let intensity = self.fragments.intensities[i] / max_intensity;
            fragments.insert((kind, self.fragments.charges[i], self.fragments.fragment_ordinals[i]), intensity);
        }

        fragments
    }

    pub fn get_observed_intensities_re_indexed(&self) -> Vec<f32> {
        // vector should be of a fixed size 174,
        let mut observed_intensities = vec![0.0; 174];
        let intensity_map = self.observed_intensity_to_fragments_map();

        for ((kind, charge, ordinal), intensity) in intensity_map.iter() {
            let ordinal_index = ordinal - 1;
            let charge_index = charge - 1;
            let charge_block = charge_index * 29;
            let kind_block = (kind * 29 * 3) as i32;
            let idx = (kind_block + charge_block + ordinal_index) as usize;
            observed_intensities[idx] = *intensity;
        }
        observed_intensities
    }

    /// Get the prosit intensities re-indexed to match the observed intensities
    ///
    /// # Arguments
    ///
    /// * `reduce_matched` - A boolean indicating whether to reduce the prosit intensities to only the matched ions
    ///
    /// # Returns
    ///
    /// * `Vec<f32>` - A vector of f32 representing the re-indexed prosit intensities, possibly reduced to only the matched ions
    pub fn get_prosit_intensities_re_indexed(&self, reduce_matched: bool) -> Vec<f32> {
        // vector should be of a fixed size 174,
        let mut prosit_intensities = vec![0.0; 174];
        let mut intensity_map = self.prosit_intensity_to_fragments_map();

        // Reduce the prosit intensities to only the matched ions
        if reduce_matched {
            let real_intensity_map = self.observed_intensity_to_fragments_map();

            for (key, _) in intensity_map.clone().iter() {
                if !real_intensity_map.contains_key(&key) {
                    intensity_map.remove(&key);
                }
            }
        }

        for ((kind, charge, ordinal), intensity) in intensity_map.iter() {
            let ordinal_index = ordinal - 1;
            let charge_index = charge - 1;
            let charge_block = charge_index * 29;
            let kind_block = (kind * 29 * 3) as i32;
            let idx = (kind_block + charge_block + ordinal_index) as usize;
            prosit_intensities[idx] = *intensity;
        }
        prosit_intensities
    }

    pub fn spectral_entropy_similarity(&self, epsilon: f32, reduce_matched: bool) -> f32 {
        let observed_intensities = self.get_observed_intensities_re_indexed();
        let prosit_intensities = self.get_prosit_intensities_re_indexed(reduce_matched);
        spectral_entropy_similarity(&observed_intensities, &prosit_intensities, epsilon)
    }

    pub fn pearson_correlation(&self, epsilon: f32, reduce_matched: bool) -> f32 {
        let observed_intensities = self.get_observed_intensities_re_indexed();
        let prosit_intensities = self.get_prosit_intensities_re_indexed(reduce_matched);
        pearson_correlation(&observed_intensities, &prosit_intensities, epsilon)
    }

    pub fn spearman_correlation(&self, epsilon: f32, reduce_matched: bool) -> f32 {
        let observed_intensities = self.get_observed_intensities_re_indexed();
        let prosit_intensities = self.get_prosit_intensities_re_indexed(reduce_matched);
        spearman_correlation(&observed_intensities, &prosit_intensities, epsilon)
    }

    pub fn cosine_similarity(&self, epsilon: f32, reduce_matched: bool) -> Option<f32> {
        let observed_intensities = self.get_observed_intensities_re_indexed();
        let prosit_intensities = self.get_prosit_intensities_re_indexed(reduce_matched);
        cosine_similarity(&observed_intensities, &prosit_intensities, epsilon)
    }

    pub fn spectral_angle_similarity(&self, epsilon: f32, reduce_matched: bool) -> f32 {
        let cosim = self.cosine_similarity(epsilon, reduce_matched).unwrap_or(0.0);
        cosim_to_spectral_angle(cosim)
    }

    pub fn prosit_intensity_to_fragments(&self) -> Fragments {
        prosit_intensities_to_fragments(self.prosit_intensity_predicted.clone())
    }

    pub fn get_feature_vector(&self, epsilon: f32, reduce_matched: bool) -> Vec<f32> {
        vec![
            self.cosine_similarity(epsilon, reduce_matched).unwrap_or(0.0),
            self.spectral_angle_similarity(epsilon, reduce_matched),
            self.pearson_correlation(epsilon, reduce_matched),
            self.spearman_correlation(epsilon, reduce_matched),
            self.spectral_entropy_similarity(epsilon, reduce_matched),
        ]
    }
}

pub fn prosit_intensities_to_fragments(
    flat_intensities: Vec<f32>,
) -> Fragments {
    let fragments = flat_prosit_array_to_fragments_map(flat_intensities);

    let mut predicted_kinds_b: Vec<Kind> = Vec::new();
    let mut predicted_kinds_y: Vec<Kind> = Vec::new();

    let mut predicted_fragment_ordinals_b = Vec::new();
    let mut predicted_fragment_ordinals_y = Vec::new();

    let mut predicted_charges_b = Vec::new();
    let mut predicted_charges_y = Vec::new();

    let mut predicted_intensities_b = Vec::new();
    let mut predicted_intensities_y = Vec::new();

    for (key, _) in fragments.iter() {

        let (kind, charge, fragment_ordinal) = key;

        let predicted_intensity = fragments.get(key).unwrap_or(&0.0);

        let kind = match kind {
            0 => Kind::B,
            1 => Kind::Y,
            _ => panic!("Invalid kind"),
        };

        if kind == Kind::B {
            predicted_kinds_b.push(kind);
            predicted_charges_b.push(*charge);
            predicted_fragment_ordinals_b.push(*fragment_ordinal);
            predicted_intensities_b.push(*predicted_intensity);
        } else {
            predicted_kinds_y.push(kind);
            predicted_charges_y.push(*charge);
            predicted_fragment_ordinals_y.push(*fragment_ordinal);
            predicted_intensities_y.push(*predicted_intensity);
        }
    }

    // invert the order of y fragments
    predicted_kinds_y.reverse();
    predicted_charges_y.reverse();
    predicted_fragment_ordinals_y.reverse();
    predicted_intensities_y.reverse();

    let fragments = Fragments {
        charges: predicted_charges_b.iter().chain(predicted_charges_y.iter()).cloned().collect(),
        kinds: predicted_kinds_b.iter().chain(predicted_kinds_y.iter()).cloned().collect(),
        fragment_ordinals: predicted_fragment_ordinals_b.iter().chain(predicted_fragment_ordinals_y.iter()).cloned().collect(),
        intensities: predicted_intensities_b.iter().chain(predicted_intensities_y.iter()).cloned().collect(),
        mz_calculated: Vec::new(),
        mz_experimental: Vec::new(),
    };

    fragments
}