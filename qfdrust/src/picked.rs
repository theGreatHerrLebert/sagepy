use std::collections::HashMap;
use crate::dataset::{Match, MatchDataset};
use crate::psm::Psm;
use rayon::prelude::*;

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

#[derive(Clone, Debug)]
struct Competition {
    forward: f32,
    forward_match: Option<Match>, // Own the Match
    reverse: f32,
    reverse_match: Option<Match>, // Own the Match
}

impl Competition {
    fn new() -> Self {
        Self {
            forward: f32::MIN,
            forward_match: None,
            reverse: f32::MIN,
            reverse_match: None,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Row {
    pub key: (String, String), // (spectrum_idx, match_idx)
    pub decoy: bool,
    pub score: f32,
    pub q_value: f64,
}

pub fn assign_q_value(
    rows: Vec<Row>,
) -> HashMap<(String, String), f64> {
    let mut q_values: HashMap<(String, String), f64> = HashMap::new();

    // Sort rows by descending score
    let mut sorted_rows = rows;
    sorted_rows.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap());

    let mut decoy_count: f64 = 1.0;
    let mut target_count: f64 = 0.0;
    let mut q_values_list: Vec<f64> = Vec::new();

    // First pass: Calculate the raw q-values
    for row in sorted_rows.iter() {
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
    for (i, row) in sorted_rows.iter_mut().enumerate().rev() {
        let q = q_values_list[i];
        if q < q_min {
            q_min = q;
        }
        row.q_value = q_min;
        q_values.insert(row.key.clone(), row.q_value);
    }

    q_values
}

pub fn tdc_picked_peptide_match(ds: &MatchDataset) -> Vec<Match> {
    use std::collections::HashMap;

    let mut peptide_map: HashMap<String, Competition> = HashMap::new();

    for matches in ds.matches.values() {
        for m in matches {
            let key = m.match_idx.clone(); // Peptide sequence

            let entry = peptide_map.entry(key).or_insert_with(Competition::new);

            if m.decoy {
                if m.score > entry.reverse {
                    entry.reverse = m.score;
                    entry.reverse_match = Some(m.clone());
                }
            } else {
                if m.score > entry.forward {
                    entry.forward = m.score;
                    entry.forward_match = Some(m.clone());
                }
            }
        }
    }

    // Collect competitions into a Vec<Competition>
    let competitions: Vec<Competition> = peptide_map.into_values().collect();

    // Prepare rows for q-value assignment
    let rows: Vec<Row> = competitions
        .into_iter()
        .flat_map(|comp| {
            let mut entries = Vec::new();
            if let Some(fwd_match) = comp.forward_match {
                entries.push(Row {
                    key: (fwd_match.spectrum_idx.clone(), fwd_match.match_idx.clone()),
                    decoy: false,
                    score: comp.forward,
                    q_value: 1.0,
                });
            }
            if let Some(rev_match) = comp.reverse_match {
                entries.push(Row {
                    key: (rev_match.spectrum_idx.clone(), rev_match.match_idx.clone()),
                    decoy: true,
                    score: comp.reverse,
                    q_value: 1.0,
                });
            }
            entries
        })
        .collect();

    // Assign q-values
    let q_values = assign_q_value(rows);

    // Update matches with q-values
    let mut result = Vec::new();
    for matches in ds.matches.values() {
        for m in matches {
            let key = (m.spectrum_idx.clone(), m.match_idx.clone());
            if let Some(&q_value) = q_values.get(&key) {
                let mut m_clone = m.clone();
                m_clone.q_value = Some(q_value);
                result.push(m_clone);
            }
        }
    }

    result
}

pub fn tdc_picked_protein_match(ds: &MatchDataset) -> Vec<Match> {
    use std::collections::HashMap;

    let mut protein_map: HashMap<String, Competition> = HashMap::new();

    for matches in ds.matches.values() {
        for m in matches {
            if let Some(protein_ids) = &m.match_identity_candidates {
                for protein_id in protein_ids {
                    let entry = protein_map.entry(protein_id.clone()).or_insert_with(Competition::new);

                    if m.decoy {
                        if m.score > entry.reverse {
                            entry.reverse = m.score;
                            entry.reverse_match = Some(m.clone());
                        }
                    } else {
                        if m.score > entry.forward {
                            entry.forward = m.score;
                            entry.forward_match = Some(m.clone());
                        }
                    }
                }
            }
        }
    }

    // Collect competitions into a Vec<Competition>
    let competitions: Vec<Competition> = protein_map.into_values().collect();

    // Prepare rows for q-value assignment
    let rows: Vec<Row> = competitions
        .into_iter()
        .flat_map(|comp| {
            let mut entries = Vec::new();
            if let Some(fwd_match) = comp.forward_match {
                entries.push(Row {
                    key: (fwd_match.spectrum_idx.clone(), fwd_match.match_idx.clone()),
                    decoy: false,
                    score: comp.forward,
                    q_value: 1.0,
                });
            }
            if let Some(rev_match) = comp.reverse_match {
                entries.push(Row {
                    key: (rev_match.spectrum_idx.clone(), rev_match.match_idx.clone()),
                    decoy: true,
                    score: comp.reverse,
                    q_value: 1.0,
                });
            }
            entries
        })
        .collect();

    // Assign q-values
    let q_values = assign_q_value(rows);

    // Update matches with q-values
    let mut result = Vec::new();
    for matches in ds.matches.values() {
        for m in matches {
            let key = (m.spectrum_idx.clone(), m.match_idx.clone());
            if let Some(&q_value) = q_values.get(&key) {
                let mut m_clone = m.clone();
                m_clone.q_value = Some(q_value);
                result.push(m_clone);
            }
        }
    }

    result
}