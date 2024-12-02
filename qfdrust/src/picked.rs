use std::collections::HashMap;
use itertools::Itertools;
use crate::psm::Psm;
use rayon::prelude::*;

#[derive(Default)]
struct Row {
    ix: String,
    decoy: bool,
    score: f32,
    q: f32,
}

#[derive(Clone, Debug)]
struct Competition {
    ix: Option<String>,
    forward: f32,
    reverse: f32,
}

impl Default for Competition {
    fn default() -> Self {
        Competition {
            ix: None,
            forward: f32::MIN,
            reverse: f32::MIN,
        }
    }
}

pub fn spectrum_q_value(scores: &Vec<Psm>, use_hyper_score: bool) -> Vec<f32> {

    // create a collection of PSMs sorted by score and keep the index
    let mut indexed_inner_collection: Vec<(usize, Psm)> = scores.iter()
        .enumerate()
        .map(|(index, item)| (index, item.clone()))
        .collect();

    // sort either by hyperscore or PSM re_score
    match use_hyper_score {
        // Sort by hyperscore
        true => {
            indexed_inner_collection.par_sort_unstable_by(|(_, a), (_, b)| b.sage_feature.hyperscore.total_cmp(&a.sage_feature.hyperscore));
        }
        // Sort by PSM re_score
        false => {
            indexed_inner_collection.par_sort_unstable_by(|(_, a), (_, b)| b.re_score.unwrap().total_cmp(&a.re_score.unwrap()));
        }
    }

    // Calculate the spectrum q-value
    let mut decoy = 1;
    let mut target = 0;

    for (_, psm) in indexed_inner_collection.iter_mut() {
        match psm.sage_feature.label == -1 {
            true => decoy += 1,
            false => target += 1,
        }
        psm.sage_feature.spectrum_q = decoy as f32 / target as f32;
    }

    // Reverse slice, and calculate the cumulative minimum
    let mut q_min = 1.0f32;
    for (_, psm) in indexed_inner_collection.iter_mut().rev() {
        q_min = q_min.min(psm.sage_feature.spectrum_q);
        psm.sage_feature.spectrum_q = q_min;
    }

    // sort the q_values by the original index
    let mut q_values = vec![0.0; scores.len()];
    for (sorted_index, psm) in indexed_inner_collection.iter() {
        q_values[*sorted_index] = psm.sage_feature.spectrum_q;
    }

    q_values
}

pub fn picked_peptide(features: &mut Vec<Psm>, use_hyper_score: bool) -> HashMap<String, f64> {

    let mut map: HashMap<String, Competition> = HashMap::default();

    for feat in features.iter() {

        let peptide_sequence_key = match feat.sage_feature.label == -1 {
            true => feat.sequence_decoy.clone().unwrap().sequence,
            false => feat.sequence.clone().unwrap().sequence,
        };

        let entry = map.entry(peptide_sequence_key).or_default();

        match feat.sage_feature.label == -1 {
            true => {
                match use_hyper_score {
                    true => {
                        entry.reverse = entry.reverse.max(feat.sage_feature.hyperscore as f32);
                    }
                    false => {
                        entry.reverse = entry.reverse.max(feat.re_score.unwrap() as f32);
                    }
                }
            }
            false => {
                match use_hyper_score {
                    true => {
                        entry.forward = entry.forward.max(feat.sage_feature.hyperscore as f32);
                    }
                    false => {
                        entry.forward = entry.forward.max(feat.re_score.unwrap() as f32);
                    }
                }
            }
        }
    }

    let q_value_map = assign_q_value(map);

    q_value_map
}

pub fn picked_protein(features: &mut Vec<Psm>, use_hyper_score: bool) -> HashMap<String, f64> {

    let mut map: HashMap<String, Competition> = HashMap::default();

    for feat in features.iter() {

        let protein_key = protein_id_from_psm(feat, "rev_", true);

        let entry = map.entry(protein_key).or_default();

        match feat.sage_feature.label == -1 {
            true => {
                match use_hyper_score {
                    true => {
                        entry.reverse = entry.reverse.max(feat.sage_feature.hyperscore as f32);
                    }
                    false => {
                        entry.reverse = entry.reverse.max(feat.re_score.unwrap() as f32);
                    }
                }
            }
            false => {
                match use_hyper_score {
                    true => {
                        entry.forward = entry.forward.max(feat.sage_feature.hyperscore as f32);
                    }
                    false => {
                        entry.forward = entry.forward.max(feat.re_score.unwrap() as f32);
                    }
                }
            }
        }
    }

    let q_value_map = assign_q_value(map);

    q_value_map
}

pub fn protein_id_from_psm(psm: &Psm, decoy_tag: &str, generate_decoys: bool) -> String {
    if psm.sage_feature.label == -1 {
        psm.proteins
            .iter()
            .map(|s| {
                if generate_decoys {
                    format!("{}{}", decoy_tag, s)
                } else {
                    s.to_string()
                }
            })
            .join(";")
    } else {
        psm.proteins.iter().join(";")
    }
}

fn assign_q_value(
    scores: HashMap<String, Competition>,
) -> HashMap<String, f64> {

    let mut q_values: HashMap<String, f64> = HashMap::new();

    let mut scores = scores
        .into_par_iter()
        .flat_map(|(_, comp)| {
            [
                (comp.ix.clone(), false, comp.forward),
                (comp.ix, true, comp.reverse),
            ]
        })
        .filter_map(|(ix, decoy, score)| {
            ix.map(|ix| Row {
                ix,
                decoy,
                score,
                q: 1.0,
            })
        })
        .collect::<Vec<Row>>();

    scores.par_sort_by(|a, b| b.score.total_cmp(&a.score));

    let mut decoy_count: f64 = 1.0;
    let mut target_count: f64 = 0.0;
    let mut q_values_list: Vec<f64> = Vec::new();

    // First pass: Calculate the raw q-values
    for row in scores.iter() {
        if row.decoy {
            decoy_count += 1.0;
        } else {
            target_count += 1.0;
        }

        // Avoid division by zero
        if target_count == 0.0 {
            q_values_list.push(1.0);
            continue;
        }

        let q = decoy_count / target_count;
        q_values_list.push(q);
    }

    // Second pass: Compute the cumulative minimum from the end
    let mut q_min = 1.0;
    for (i, row) in scores.iter_mut().enumerate().rev() {
        let q = q_values_list[i];
        if q < q_min {
            q_min = q;
        }
        row.q = q_min as f32;
        q_values.insert(row.ix.clone(), row.q as f64);
    }

    q_values
}