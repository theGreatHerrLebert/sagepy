use rand::prelude::*;

/// Use target-decoy competition to calculate q-values
///
/// # Arguments
///
/// * `scores` - A vector of floats representing the scores
/// * `target` - A vector of booleans representing the target/decoy status
/// * `desc` - A boolean representing the sort order of the scores
///
/// # Returns
///
/// * `Vec<f64>` - A vector of floats representing the q-values
///
pub fn target_decoy_competition(scores: &Vec<f64>, target: &Vec<bool>, desc: bool) -> Vec<f64> {
    assert_eq!(scores.len(), target.len(), "Scores and target must be the same length");

    // Create a vector of indices and sort by scores
    let mut indices: Vec<usize> = (0..scores.len()).collect();
    if desc {
        indices.sort_by(|&i, &j| scores[j].partial_cmp(&scores[i]).unwrap());
    } else {
        indices.sort_by(|&i, &j| scores[i].partial_cmp(&scores[j]).unwrap());
    }

    // Apply sorted indices to scores and targets
    let sorted_scores: Vec<f64> = indices.iter().map(|&i| scores[i]).collect();
    let sorted_target: Vec<bool> = indices.iter().map(|&i| target[i]).collect();

    // Calculate cumulative sums for targets and decoys
    let mut cum_targets = 0;
    let mut cum_decoys = 0;
    let mut cum_targets_vec = Vec::new();
    let mut cum_decoys_vec = Vec::new();

    for &t in &sorted_target {
        if t {
            cum_targets += 1;
        } else {
            cum_decoys += 1;
        }
        cum_targets_vec.push(cum_targets);
        cum_decoys_vec.push(cum_decoys);
    }

    // Calculate FDR
    let mut fdr: Vec<f64> = cum_decoys_vec.iter()
        .zip(cum_targets_vec.iter())
        .map(|(&d, &t)| if t > 0 { (d as f64 + 1.0) / t as f64 } else { 1.0 })
        .collect();

    // Calculate q-values
    fdr.reverse();
    let reversed_scores: Vec<f64> = sorted_scores.iter().rev().cloned().collect();
    let mut q_vals = fdr_to_q_value(&reversed_scores, &fdr);
    q_vals.reverse();

    // Reorder q_vals to original order
    let mut final_q_vals = vec![0.0; scores.len()];
    for (original_pos, &sorted_pos) in indices.iter().enumerate() {
        final_q_vals[sorted_pos] = q_vals[original_pos];
    }

    final_q_vals
}

/// Convert FDR to q-values
///
/// # Arguments
///
/// * `scores` - A vector of floats representing the scores
/// * `fdr` - A vector of floats representing the FDR
///
/// # Returns
///
/// * `Vec<f64>` - A vector of floats representing the q-values
///
fn fdr_to_q_value(scores: &[f64], fdr: &[f64]) -> Vec<f64> {
    assert_eq!(scores.len(), fdr.len(), "Scores and FDR must be of the same length");

    let mut min_q = 1.0;
    let mut qvals = vec![1.0; fdr.len()];
    let mut start = 0;

    for (idx, &score) in scores.iter().enumerate() {
        // check if the next score is the same
        if idx < scores.len() - 1 && scores[idx + 1] == score {
            continue;
        }

        // update the minimum q-value
        if fdr[start] < min_q {
            min_q = fdr[start];
        }

        for qval in &mut qvals[start..=idx] {
            *qval = min_q;
        }
        start = idx + 1;
    }

    qvals
}

fn _estimate_pi0(pval_list: &Vec<f64>) -> f64 {
    let num_lambda = 100;
    let max_lambda = 0.5;
    let num_boot = 100;
    let max_size = 1000;
    let mut rng = rand::rng();

    let n_pval = pval_list.len();
    let mut pi0s_list = Vec::new();
    let mut lambda_list = Vec::new();

    for idx in 0..num_lambda {
        let cur_lambda = ((idx + 1) as f64 / num_lambda as f64) * max_lambda;
        let start = pval_list.binary_search_by(|p| p.partial_cmp(&cur_lambda).unwrap()).unwrap_or_else(|pos| pos);
        let w1 = n_pval - start;
        let pi0 = w1 as f64 / n_pval as f64 / (1.0 - cur_lambda);

        if pi0 > 0.0 {
            lambda_list.push(cur_lambda);
            pi0s_list.push(pi0);
        }
    }

    assert!(!pi0s_list.is_empty(), "Error in the input data: too good separation between target and decoy PSMs.");

    let min_pi0 = *pi0s_list.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
    let mut mse_list = vec![0.0; pi0s_list.len()];

    let dist = rand::distr::Uniform::new(0, n_pval).unwrap();

    for _ in 0..num_boot {
        let num_draw = std::cmp::min(n_pval, max_size);
        let mut p_boot_list: Vec<f64> = (0..num_draw).map(|_| pval_list[dist.sample(&mut rng)]).collect();
        p_boot_list.sort_by(|a, b| a.partial_cmp(b).unwrap());

        for (idx, &lambda) in lambda_list.iter().enumerate() {
            let start = p_boot_list.binary_search_by(|p| p.partial_cmp(&lambda).unwrap()).unwrap_or_else(|pos| pos);
            let w1 = num_draw - start;
            let pi0_boot = w1 as f64 / num_draw as f64 / (1.0 - lambda);
            mse_list[idx] += (pi0_boot - min_pi0).powi(2);
        }
    }

    let min_idx = mse_list.iter().enumerate().min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap()).map(|(idx, _)| idx).unwrap();

    pi0s_list[min_idx].clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn setup_desc_scores() -> (Vec<f64>, Vec<bool>, Vec<f64>) {
        let scores = vec![10.0, 10.0, 9.0, 8.0, 7.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 2.0, 1.0, 1.0, 1.0, 1.0];
        let target = vec![true, true, true, true, false, true, true, false, true, false, true, false, false, false, false, false];
        let q_vals = vec![0.25, 0.25, 0.25, 0.25, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.42857142857142855, 0.42857142857142855, 0.5714285714285714, 0.625, 0.625, 1.0, 1.0, 1.0, 1.0];
        (scores, target, q_vals)
    }

    #[test]
    fn test_tdc_descending() {
        let (scores, target, true_q_vals) = setup_desc_scores();
        let q_vals = target_decoy_competition(&scores, &target, true);
        assert_eq!(q_vals, true_q_vals, "Q-values for descending scores are incorrect.");
    }

    #[test]
    fn test_tdc_ascending() {
        let (mut scores, target, true_q_vals) = setup_desc_scores();
        scores = scores.into_iter().map(|x| -x).collect(); // Negate scores for ascending test
        let q_vals = target_decoy_competition(&scores, &target, false);
        assert_eq!(q_vals, true_q_vals, "Q-values for ascending scores are incorrect.");
    }
}