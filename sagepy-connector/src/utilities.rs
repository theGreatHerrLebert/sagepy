use unimod::unimod::{quanzie_mass, quantized_mass_to_unimod};
use std::collections::HashSet;

/// Convert a Sage sequence and modifications to a Unimod sequence
///
/// # Arguments
///
/// * `sequence` - A string representing the amino acid sequence
/// * `modifications` - A vector of floats representing the modifications
///
/// # Returns
///
/// * `String` - A string representing the Unimod sequence
///
pub fn sage_sequence_to_unimod_sequence(sequence: String, modifications: &Vec<f32>, expected_modifications: &HashSet<String>) -> String {

    assert_eq!(sequence.len(), modifications.len(), "Sequence and modifications must be the same length");

    // go over each char and check if modification is present (not 0.0) and possibly convert to unimod
    let mut unimod_sequence = String::new();
    let unimod_modifications_qunatized = quantized_mass_to_unimod();
    let empty_vec = Vec::new();

    for (idx, aa) in sequence.chars().enumerate() {

        // add amino acid to the unimod sequence
        unimod_sequence.push(aa);

        // check if the modification is nonzero, need to translate to unimod
        if modifications[idx] != 0.0 {

            // quantize the mass from nonzero modification
            let quantized_mass = quanzie_mass(modifications[idx]);

            // find the candidate modifications for the quantized mass
            let modifications = unimod_modifications_qunatized.get(&quantized_mass).unwrap_or(&empty_vec);

            let mut found = false;

            // check if the expected modification is in the candidate modifications
            for modification in modifications {
                if expected_modifications.contains(&modification.to_string()) {
                    unimod_sequence.push_str(modification);
                    found = true;
                }
            }

            // if the expected modification is not found, add a placeholder
            if !found {
                unimod_sequence.push_str("[UNIMOD:X]");
            }
        }
    }
    unimod_sequence
}