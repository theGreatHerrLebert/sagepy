from typing import List, Union, Dict
from sagepy.core.scoring import Psm
from sagepy.core import IndexedDatabase

import sagepy_connector
psc = sagepy_connector.py_mobility_model

def predict_sage_im(
        psm_collection: Union[List[Psm], Dict[str, List[Psm]]],
        indexed_db: IndexedDatabase) -> None:
    """ Predict ion mobility using SAGE IM model.
    Args:
        psm_collection: a list of features
        indexed_db: an indexed database
    """

    f_collection = []

    if isinstance(psm_collection, dict):
        for _, values in psm_collection.items():
            f_collection.extend(values)

    else:
        f_collection = psm_collection

    psc.py_predict_im(
        [feature.get_py_ptr() for feature in f_collection],
        indexed_db.get_py_ptr(),
    )