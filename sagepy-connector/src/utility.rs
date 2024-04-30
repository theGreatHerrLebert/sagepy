use sage_core::ion_series::Kind;
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

pub fn flat_prosit_array_to_py_fragments(flat_intensities: Vec<f32>) -> PyFragments {

    let reshaped_intensities = reshape_prosit_array(flat_intensities);
    let mut kinds: Vec<Kind> = Vec::new();
    let mut ordinals: Vec<u32> = Vec::new();
    let mut intensities: Vec<Vec<f32>> = Vec::new();

    for z in 1..=3 {
        let intensity_b: Vec<f32> = reshaped_intensities[..].iter().map(|x| x[1][z as usize - 1]).collect();
        let intensity_y: Vec<f32> = reshaped_intensities[..].iter().map(|x| x[0][z as usize - 1]).collect(); // Reverse for y
    }

    todo!("Finish this function")
}