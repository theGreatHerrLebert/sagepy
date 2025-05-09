use pyo3::prelude::*;
use pyo3::types::PyList;
use pyo3::exceptions::PyRuntimeError;
use sage_core::ml::retention_model::predict;
use sage_core::scoring::Feature;
use crate::py_database::PyIndexedDatabase;
use crate::py_scoring::{PyPsm};

#[pyfunction]
pub fn py_predict_rt(
    _py: Python,
    psm_collection: &Bound<'_, PyList>,
    indexed_database: &PyIndexedDatabase,
) -> PyResult<()> {
    
    let indexed_feats: Vec<(usize, Feature)> = psm_collection.iter()
        .enumerate()
        .map(|(idx, item)| {
            let psm: Bound<'_, PyPsm> = item
                .extract()
                .expect("Failed to extract PyPsm");
            // clone just the inner Feature (sage_feature)
            (idx, psm.borrow().inner.sage_feature.clone())
        })
        .collect();
    
    let mut feats: Vec<Feature> = indexed_feats.iter()
        .map(|(_, feat)| feat.clone())
        .collect();

    if predict(&indexed_database.inner, &mut feats).is_none() {
        return Err(PyRuntimeError::new_err(
            "Retention model fit failed: not enough data or R² < 0.7"
        ));
    }

    // 3) write back the two mutated fields
    for ((orig_idx, _), updated) in indexed_feats.iter().zip(feats.iter()) {
        let psm: Bound<'_, PyPsm> = psm_collection
            .get_item(*orig_idx)
            .expect("Failed to get PyPsm")
            .extract()?;
        let mut psm_borrow = psm.borrow_mut();
        psm_borrow.inner.sage_feature.predicted_rt   = updated.predicted_rt;
        psm_borrow.inner.sage_feature.delta_rt_model = updated.delta_rt_model;
    }

    Ok(())
}


#[pymodule]
pub fn py_retention_time_prediction(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(py_predict_rt, m)?)?;
    Ok(())
}