use unimod::unimod::{quanzie_mass, quantized_mass_to_unimod};
use std::collections::HashSet;

/// Convert a Sage sequence and modifications to a Unimod/ProForma-like sequence
///
/// # Arguments
///
/// * `sequence` - A string representing the amino acid sequence
/// * `modifications` - A vector of floats representing the modifications (per-residue delta mass)
///
/// # Returns
///
/// * `String` - A string representing the sequence with UNIMOD or delta-mass annotations
///
pub fn sage_sequence_to_unimod_sequence(
    sequence: String,
    modifications: &Vec<f32>,
    expected_modifications: &HashSet<String>,
) -> String {
    assert_eq!(
        sequence.len(),
        modifications.len(),
        "Sequence and modifications must be the same length"
    );

    let mut unimod_sequence = String::new();
    let unimod_modifications_qunatized = quantized_mass_to_unimod();
    let empty_vec: Vec<&str> = Vec::new();

    for (idx, aa) in sequence.chars().enumerate() {
        unimod_sequence.push(aa);

        let mass_shift = modifications[idx];
        if mass_shift != 0.0 {
            // quantize mass to lookup bucket
            let quantized_mass = quanzie_mass(mass_shift);

            // candidate UNIMOD strings (likely like "[UNIMOD:5]")
            let candidates = unimod_modifications_qunatized
                .get(&quantized_mass)
                .unwrap_or(&empty_vec);

            // collect all candidates that are in expected_modifications, stripping outer brackets
            let mut matching: Vec<String> = candidates
                .iter()
                .filter(|m| expected_modifications.contains(&m.to_string()))
                .map(|m| {
                    let s = m.to_string();
                    s.trim_start_matches('[')
                        .trim_end_matches(']')
                        .to_string()
                })
                .collect();

            if !matching.is_empty() {
                // stable output
                matching.sort_unstable();

                unimod_sequence.push('[');
                unimod_sequence.push_str(&matching.join("|"));
                unimod_sequence.push(']');
            } else {
                // unknown → print delta mass with fixed precision
                unimod_sequence.push('[');
                unimod_sequence.push_str(&format!("{:+.5}", mass_shift));
                unimod_sequence.push(']');
            }
        }
    }

    unimod_sequence
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_unmodified_sequence() {
        let sequence = "PEPTIDE".to_string();
        let modifications = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let expected_mods = HashSet::new();

        let result = sage_sequence_to_unimod_sequence(sequence, &modifications, &expected_mods);
        assert_eq!(result, "PEPTIDE");
    }

    #[test]
    fn test_unknown_modification_shows_delta_mass() {
        let sequence = "PEP".to_string();
        let modifications = vec![0.0, 15.9949, 0.0]; // Oxidation mass
        let expected_mods = HashSet::new(); // Empty - no expected mods

        let result = sage_sequence_to_unimod_sequence(sequence, &modifications, &expected_mods);
        // Should show delta mass since no expected mods match
        assert!(result.starts_with("PE[+"));
        assert!(result.contains("15.99"));
    }

    #[test]
    fn test_known_modification_shows_unimod() {
        let sequence = "PEP".to_string();
        let modifications = vec![0.0, 15.9949, 0.0]; // Oxidation mass
        let mut expected_mods = HashSet::new();
        expected_mods.insert("[UNIMOD:35]".to_string()); // Oxidation

        let result = sage_sequence_to_unimod_sequence(sequence, &modifications, &expected_mods);
        assert_eq!(result, "PE[UNIMOD:35]P");
    }

    #[test]
    #[should_panic(expected = "Sequence and modifications must be the same length")]
    fn test_mismatched_lengths_panics() {
        let sequence = "PEPTIDE".to_string();
        let modifications = vec![0.0, 0.0]; // Too short
        let expected_mods = HashSet::new();

        sage_sequence_to_unimod_sequence(sequence, &modifications, &expected_mods);
    }
}