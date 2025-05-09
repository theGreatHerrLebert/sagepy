import os
from typing import List, Union, Dict, Any
import pandas as pd

from imspy.timstof.dda import TimsDatasetDDA

from imspy.timstof.dbsearch.utility import (
    get_ms1_ims_spectrum,
    sanitize_mz,
    sanitize_charge,
    get_searchable_spec
)

from sagepy.core import (
    SpectrumProcessor,
    Precursor,
    Tolerance
)

def process_timstof_datasets(
    dataset_dirs: Union[str, List[str]],
    use_bruker_sdk: bool = False,
    max_peaks: int = 10_000,
    num_threads: int = 16,
    ms1_take_top_n: int = 10_000,
    ms1_deisotope: bool = True,
    fragment_take_top_n: int = 150
) -> Dict[str, Dict[str, Any]]:
    """
    Process one or more Bruker TIMS .d folders to extract summarized fragment DataFrames
    and MS1 spectra suitable for downstream search.

    Parameters
    ----------
    dataset_dirs : str or list of str
        A path or list of paths to the .d folders to process.
    use_bruker_sdk : bool, default False
        Whether to use the Bruker SDK when reading data.
    max_peaks : int, default 100000
        Maximum number of peaks to collect per precursor frame.
    num_threads : int, default 16
        Number of threads for parallel data extraction.
    ms1_take_top_n : int, default 10000
        Number of top peaks to keep when processing MS1 spectra.
    ms1_deisotope : bool, default True
        Whether to deisotope MS1 spectra.
    fragment_take_top_n : int, default 150
        Number of top peaks to keep when processing fragment spectra.

    Returns
    -------
    results : dict
        A mapping from each dataset directory to a dict with keys:
            'fragments'   : pd.DataFrame of summarized fragment ions,
            'ms1_spectra' : List of processed MS1 spectra objects.
    """
    if isinstance(dataset_dirs, str):
        dataset_dirs = [dataset_dirs]

    results: Dict[str, Dict[str, Any]] = {}
    for file_id, dataset_dir in enumerate(dataset_dirs):
        ds_name = os.path.basename(dataset_dir.rstrip(os.sep))
        handle = TimsDatasetDDA(dataset_dir, use_bruker_sdk=use_bruker_sdk)

        # Extract precursor frames and process MS1 spectra
        precursor_frames = handle.get_precursor_frames(
            max_peaks=max_peaks,
            num_threads=num_threads
        )
        ms1_processor = SpectrumProcessor(
            take_top_n=ms1_take_top_n,
        )
        ms1_spectra = [
            get_ms1_ims_spectrum(
                raw_spectrum=spec,
                spec_id=f"{spec.frame_id}-{ds_name}",
                time=spec.retention_time / 60,
                spec_processor=ms1_processor,
                file_id=file_id
            ) for spec in precursor_frames
        ]

        # Extract and summarize PASEF fragments
        fragments = handle.get_pasef_fragments(num_threads=num_threads)
        fragments = fragments.groupby('precursor_id').agg({
            'frame_id': 'first',
            'time': 'first',
            'precursor_id': 'first',
            'raw_data': 'sum',
            'scan_begin': 'first',
            'scan_end': 'first',
            'isolation_mz': 'first',
            'isolation_width': 'first',
            'collision_energy': 'first',
            'largest_peak_mz': 'first',
            'average_mz': 'first',
            'monoisotopic_mz': 'first',
            'charge': 'first',
            'average_scan': 'first',
            'intensity': 'first',
            'parent_id': 'first',
        })

        # Compute marginal ion mobility and build spec IDs
        fragments['mobility'] = fragments['raw_data'].apply(
            lambda rd: rd.get_inverse_mobility_along_scan_marginal()
        )
        fragments['spec_id'] = fragments.apply(
            lambda r: f"{r.frame_id}-{r.precursor_id}-{ds_name}",
            axis=1
        )

        # Build Precursor objects
        fragments['sage_precursor'] = fragments.apply(
            lambda r: Precursor(
                mz=sanitize_mz(r['monoisotopic_mz'], r['largest_peak_mz']),
                intensity=r['intensity'],
                charge=sanitize_charge(r['charge']),
                isolation_window=Tolerance(da=(-3, 3)),
                collision_energy=r['collision_energy'],
                inverse_ion_mobility=r['mobility'],
                spectrum_ref=r['spec_id']
            ),
            axis=1
        )

        # Process fragment spectra for searching
        fragments['processed_spec'] = fragments.apply(
            lambda r: get_searchable_spec(
                precursor=r['sage_precursor'],
                raw_fragment_data=r['raw_data'],
                spec_processor=SpectrumProcessor(take_top_n=fragment_take_top_n),
                spec_id=r['spec_id'],
                time=r['time'],
                file_id=file_id
            ),
            axis=1
        )

        results[dataset_dir] = {
            'fragments': fragments,
            'ms1_spectra': ms1_spectra
        }

    return results