use std::collections::HashMap;
use sage_core::ion_series::Kind;
use sage_core::scoring::Fragments;
use crate::py_scoring::PyFragments;

/// Calculates the cosine similarity between two vectors.
///
/// # Arguments
///
/// * `vec1` - A vector of f32.
/// * `vec2` - A vector of f32.
///
/// # Returns
///
/// * The cosine similarity as a f32. If vectors are empty or different lengths, it returns None.
///
pub fn cosine_similarity(vec1: &Vec<f32>, vec2: &Vec<f32>) -> Option<f32> {
    if vec1.len() != vec2.len() || vec1.is_empty() {
        return None;
    }

    let dot_product: f32 = vec1.iter().zip(vec2.iter()).map(|(a, b)| a * b).sum();
    let magnitude_vec1: f32 = vec1.iter().map(|x| x.powi(2)).sum::<f32>().sqrt();
    let magnitude_vec2: f32 = vec2.iter().map(|x| x.powi(2)).sum::<f32>().sqrt();

    if magnitude_vec1 == 0.0 || magnitude_vec2 == 0.0 {
        return None;
    }

    Some(dot_product / (magnitude_vec1 * magnitude_vec2))
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

pub fn flat_prosit_array_to_py_fragments(flat_intensities: Vec<f32>) -> HashMap<(u32, i32, i32), f32> {
    // Reshape the flat prosit array into a 3D array of shape (29, 2, 3)
    let reshaped_intensities = reshape_prosit_array(flat_intensities);

    // create hashmap of (kind, ordinal, charge) -> intensity
    let mut fragments: HashMap<(u32, i32, i32), f32> = HashMap::new();
    for z in 1..=3 {
        let intensity_b: Vec<f32> = reshaped_intensities[..].iter().map(|x| x[1][z as usize - 1]).collect();
        for i in 1..=29 {
                let intensity = intensity_b[i as usize - 1];
                if intensity >= 0.0 {
                    fragments.insert((0, i, z), intensity);
                }
        }

        let mut intensity_y: Vec<f32> = reshaped_intensities[..].iter().map(|x| x[0][z as usize - 1]).collect();
        intensity_y.reverse();
        for i in 1..=29 {
            let intensity = intensity_y[i as usize - 1];
            if intensity >= 0.0 {
                fragments.insert((1, i, z), intensity);
            }
        }
    }
    fragments
}

pub fn _map_to_py_fragments(fragments: &HashMap<(u32, i32, i32), f32>,
                           mz_calculated: Vec<f32>, mz_experimental: Vec<f32>) -> PyFragments {

    let mut kinds: Vec<Kind> = Vec::new();
    let mut ordinals: Vec<i32> = Vec::new();
    let mut charges: Vec<i32> = Vec::new();
    let mut intensities: Vec<f32> = Vec::new();

    for (kind, ordinal, charge) in fragments.keys() {
        let intensity = fragments.get(&(*kind, *ordinal, *charge)).unwrap();
        let kind = match kind {
            0 => Kind::B,
            1 => Kind::Y,
            _ => panic!("Invalid ion kind"),
        };
        kinds.push(kind);
        ordinals.push(*ordinal);
        charges.push(*charge);
        intensities.push(*intensity);
    }

    let fragments = Fragments {
        mz_calculated,
        mz_experimental,
        kinds,
        fragment_ordinals: ordinals,
        charges,
        intensities,
    };

    PyFragments {
        inner: fragments,
    }
}

pub fn py_fragments_to_map(fragments: &PyFragments) -> HashMap<(u32, i32, i32), f32> {
    let mut fragments_map: HashMap<(u32, i32, i32), f32> = HashMap::new();
    for i in 0..fragments.inner.mz_calculated.len() {
        let kind = match fragments.inner.kinds[i] {
            Kind::B => 0,
            Kind::Y => 1,
            _ => panic!("Invalid ion kind"),
        };
        fragments_map.insert((kind, fragments.inner.fragment_ordinals[i],
                              fragments.inner.charges[i]), fragments.inner.intensities[i]);
    }
    fragments_map
}