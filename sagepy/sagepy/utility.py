from typing import Dict, Tuple, List, Optional, Union

import re

from numba import jit
import numpy as np
import pandas as pd
from numpy.typing import NDArray

from sagepy.core.scoring import Psm
from sagepy.core.spectrum import ProcessedSpectrum, RawSpectrum, Precursor, SpectrumProcessor, Representation
from sagepy.core.mass import Tolerance
from sagepy.core.database import IndexedDatabase, EnzymeBuilder, SageSearchConfiguration
from sagepy.qfdr.tdc import target_decoy_competition_pandas

from pyteomics import mzml

import sagepy_connector
psc = sagepy_connector.py_utility
from typing import Iterator


@jit(nopython=True)
def calculate_ppm_error(measured_value, reference_value):
    ppm_error = ((measured_value - reference_value) / reference_value) * 1_000_000
    return ppm_error


@jit(nopython=True)
def calculate_ppms(measured_values, reference_values):
    n = len(measured_values)
    ppm_errors = np.empty(n, dtype=np.float64)
    for i in range(n):
        ppm_errors[i] = calculate_ppm_error(measured_values[i], reference_values[i])
    return ppm_errors


@jit(nopython=True)
def mean_ppm(mz_observed, mz_calculated) -> float:
    return np.mean(calculate_ppms(mz_observed, mz_calculated))


@jit(nopython=True)
def median_ppm(mz_observed, mz_calculated) -> float:
    return np.median(calculate_ppms(mz_observed, mz_calculated))


def prosit_intensities_to_fragments_map(intensities: NDArray) -> Dict[Tuple[int, int, int], float]:
    """ Convert a Prosit intensity array to a fragment map (ion_type, charge, ordinal) -> intensity.

    Args:
        intensities: a Prosit intensity array

    Returns:
        a fragment map (ion_type, charge, ordinal) -> intensity
    """
    fragment_map = psc.flat_prosit_array_to_fragments_map(intensities)
    return fragment_map


def py_fragments_to_fragments_map(fragments, normalize: bool = True) -> Dict[Tuple[int, int, int], float]:
    """ Convert a Fragments object to a fragment map (ion_type, charge, ordinal) -> intensity.

    Args:
        normalize: whether to normalize the intensities
        fragments: a Fragments object

    Returns:
        a fragment map (ion_type, charge, ordinal) -> intensity
    """
    return psc.py_fragments_to_fragments_map(fragments.get_py_ptr(), normalize)


def get_features(ds: pd.DataFrame, score: Optional[str] = None) -> (NDArray, NDArray):

    score = score if score is not None else "hyperscore"

    features = [
        f"{score}",
        "delta_rt",
        "delta_ims",
        "cosine_similarity",
        "delta_mass",
        "rank",
        "isotope_error",
        "average_ppm",
        "delta_next",
        "delta_best",
        "matched_peaks",
        "longest_b",
        "longest_y",
        "longest_y_pct",
        "missed_cleavages",
        "matched_intensity_pct",
        "poisson",
        "charge",
        "intensity_ms1",
        "intensity_ms2",
        "collision_energy",
        "cosine_similarity",
        "spectral_angle_similarity",
        "pearson_correlation",
        "spearman_correlation",
        "spectral_entropy_similarity",
    ]
    ds = ds.copy()

    ds["intensity_ms1"] = ds["intensity_ms1"].apply(lambda x: np.log(x + 1))
    ds["intensity_ms2"] = ds["intensity_ms2"].apply(lambda x: np.log(x + 1))

    X = ds[features].to_numpy().astype(np.float32)

    # make sure that there are no NaN values
    X = np.nan_to_num(X, nan=0.0)

    Y = ds["decoy"].to_numpy()
    Y = np.array([0 if x else 1 for x in Y]).astype(np.float32)

    return X, Y


def apply_mz_calibration(psm, fragments: pd.DataFrame, use_median: bool = True,
                         target_q: float = 0.01) -> float:
    """Apply a mass calibration to a set of peptide spectrum matches

    Args:
        psm: The peptide spectrum matches
        fragments: The fragments
        use_median: Whether to use the median ppm error
        target_q: The target q-value

    Returns:
        float: The ppm error
    """

    psms = []
    for _, item in psm.items():
        psms.extend(item)

    P = psm_collection_to_pandas(psms)
    TDC = target_decoy_competition_pandas(P, method="psm", score="hyperscore")
    TDC = TDC[TDC.q_value <= target_q]

    B = pd.merge(P, TDC, on=["spec_idx", "match_idx"])

    ppm_error = 0.0

    if use_median:
        ppm_error = np.median(B.delta_mass)
    else:
        ppm_error = np.mean(B.delta_mass)

    fragments.apply(lambda row: row.processed_spec.calibrate_mz_ppm(ppm_error), axis=1)

    return ppm_error


def create_query(
        precursor_mz: float,
        precursor_charge: Optional[int],
        precursor_intensity: float,
        isolation_window_lower: float,
        isolation_window_upper: float,
        collision_energy: float,
        retention_time: float,
        ion_injection_time: float,
        total_ion_current: float,
        fragment_mz: NDArray,
        fragment_intensity: NDArray,
        spec_id: str,
        isolation_window_in_dalton: bool = True,
        file_id: int = 0,
        ms_level: int = 2,
        take_top_n_peaks: int = 150,
        min_fragment_mz: float = 100,
        max_fragment_mz: float = 2000,
) -> ProcessedSpectrum:
    """Create a query spectrum

    Args:
        precursor_mz: The precursor m/z
        precursor_charge: The precursor charge
        precursor_intensity: The precursor intensity
        isolation_window_lower: The lower isolation window
        isolation_window_upper: The upper isolation window
        collision_energy: The collision energy
        retention_time: The retention time
        ion_injection_time: The ion injection time
        total_ion_current: The total ion current
        fragment_mz: The fragment m/z
        fragment_intensity: The fragment intensity
        spec_id: The spectrum ID
        isolation_window_in_dalton: Whether the isolation window is in Da
        file_id: The file ID
        ms_level: The MS level
        take_top_n_peaks: The number of top peaks to take
        min_fragment_mz: The minimum fragment m/z
        max_fragment_mz: The maximum fragment m/z

    Returns:
        ProcessedSpectrum: The processed spectrum
    """

    # configure the spectrum processor
    spec_processor = SpectrumProcessor(take_top_n_peaks, min_fragment_mz, max_fragment_mz)

    # set selection window bounds
    if isolation_window_in_dalton:
        tolerance = Tolerance(da=(isolation_window_lower, isolation_window_upper))
    else:
        tolerance = Tolerance(ppm=(isolation_window_lower, isolation_window_upper))

    # create the precursor that was fragmented
    sage_precursor = Precursor(
        mz=precursor_mz,
        intensity=precursor_intensity,
        charge=precursor_charge,
        isolation_window=tolerance,
        collision_energy=collision_energy,
    )

    # create the raw spectrum of fragment ions with precursor information
    spec = RawSpectrum(
        file_id=file_id,
        ms_level=ms_level,
        spec_id=spec_id,
        representation=Representation(),
        precursors=[sage_precursor],
        scan_start_time=retention_time,
        ion_injection_time=ion_injection_time,
        total_ion_current=total_ion_current,
        mz=fragment_mz.astype(np.float32),
        intensity=fragment_intensity.astype(np.float32)
    )

    # process the spectrum
    processed_spec = spec_processor.process(spec)

    return processed_spec

def create_sage_database(
    fasta_path: str,
    missed_cleavages: int = 2,
    min_len: int = 7,
    max_len: int = 30,
    cleave_at: str = 'KR',
    restrict: str = 'P',
    c_terminal: bool = True,
    generate_decoys: bool = True,
    bucket_size: int = 16384,
    static_mods: Union[Dict[str, str], Dict[int, str]] = {"C": "[UNIMOD:4]"},
    variable_mods: Union[Dict[str, List[int]], Dict[int, List[str]]] = {"M": ["[UNIMOD:35]"], "[": ["[UNIMOD:1]"]},
) -> IndexedDatabase:
    """Create a SAGE database

    Args:
        fasta_path: The path to the FASTA file
        missed_cleavages: The number of missed cleavages
        min_len: The minimum peptide length
        max_len: The maximum peptide length
        cleave_at: The cleavage sites
        restrict: The restriction sites
        c_terminal: Whether the enzyme is C-terminal
        generate_decoys: Whether to generate decoys
        bucket_size: The bucket size
        static_mods: The static modifications
        variable_mods: The variable modifications

    Returns:
        IndexedDatabase: The indexed database ready for searching
    """

    # Configure enzyme digestion
    enzyme_builder = EnzymeBuilder(
        missed_cleavages=missed_cleavages,
        min_len=min_len,
        max_len=max_len,
        cleave_at=cleave_at,
        restrict=restrict,
        c_terminal=c_terminal,
    )

    # Read FASTA file
    with open(fasta_path, 'r') as infile:
        fasta = infile.read()

    # Set up SAGE configuration
    sage_config = SageSearchConfiguration(
        fasta=fasta,
        static_mods=static_mods,
        variable_mods=variable_mods,
        enzyme_builder=enzyme_builder,
        generate_decoys=generate_decoys,
        bucket_size=bucket_size,
    )

    # Generate and return the indexed database
    indexed_db = sage_config.generate_indexed_database()
    return indexed_db


def extract_mzml_data(file_path: str) -> pd.DataFrame:
    """
    Extract relevant data from an mzML file
    Args:
        file_path: Path to the mzML file

    Returns:
        pd.DataFrame: A pandas DataFrame with the extracted data

    """
    d = []

    with mzml.read(file_path) as reader:
        for i, spectrum in enumerate(reader):
            # Check if the spectrum is an MS2 (DDA data usually has MS2 spectra)
            if spectrum['ms level'] == 2:
                precursor = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]

                spec_id = spectrum['id']
                total_ion_current = spectrum['total ion current']
                retention_time = spectrum['scanList']['scan'][0]['scan start time']
                injection_time = spectrum['scanList']['scan'][0]['ion injection time']
                collision_energy = spectrum['precursorList']['precursor'][0]['activation']['collision energy']

                # Extract the relevant metadata
                isolation_window = spectrum['precursorList']['precursor'][0]['isolationWindow']

                lower = isolation_window['isolation window target m/z'] - isolation_window[
                    'isolation window lower offset']
                upper = isolation_window['isolation window target m/z'] + isolation_window[
                    'isolation window upper offset']

                precursor_mz = precursor['selected ion m/z']
                precursor_charge = precursor['charge state']
                precursor_intensity = precursor.get('peak intensity', None)

                # Extract the fragment ion spectrum (m/z and intensity arrays)
                fragment_mz = spectrum['m/z array']
                fragment_intensity = spectrum['intensity array']

                processed_spec = create_query(
                    precursor_mz=precursor_mz,
                    precursor_charge=precursor_charge,
                    precursor_intensity=precursor_intensity,
                    isolation_window_lower=lower,
                    isolation_window_upper=upper,
                    collision_energy=collision_energy,
                    retention_time=retention_time,
                    ion_injection_time=injection_time,
                    total_ion_current=total_ion_current,
                    fragment_mz=fragment_mz,
                    fragment_intensity=fragment_intensity,
                    spec_id=spec_id,
                )

                # Create a row dictionary with relevant data
                row = {
                    "spec_id": spec_id,
                    "precursor_mz": precursor_mz,
                    "precursor_charge": precursor_charge,
                    "precursor_intensity": precursor_intensity,
                    "retention_time": retention_time,
                    "injection_time": injection_time,
                    "collision_energy": collision_energy,
                    "total_ion_current": total_ion_current,
                    "processed_spec": processed_spec
                }

                # Append the row to the list
                d.append(row)

    # Convert the list of dictionaries to a pandas DataFrame
    exp_data = pd.DataFrame(d)
    return exp_data

def psm_collection_to_feature_matrix(psm_collection: Union[List[Psm], Dict[str, List[Psm]]], num_threads: int = 4) -> NDArray:
    """Convert a list of peptide spectrum matches to a dictionary

    Args:
        psm_collection (Union[List[Psm], Dict[str, List[Psm]]): The peptide spectrum matches
        num_threads (int, optional): The number of threads to use. Defaults to 4.

    Returns:
        Dict[str, List[Psm]]: The dictionary of peptide spectrum matches
    """

    psms = []

    if isinstance(psm_collection, dict):
        for _, psm_candidates in psm_collection.items():
            psms.extend(psm_candidates)

    else:
        psms = psm_collection

    return np.array(psc.psms_to_feature_matrix([psm.get_py_ptr() for psm in psms], num_threads))

# get_psm_sequences_par

def get_psm_sequences(psm_collection: Union[List[Psm], Dict[str, List[Psm]]], num_threads: int = 4) -> List[str]:
    """Get the peptide sequences from a list of peptide spectrum matches

    Args:
        psm_collection (Union[List[Psm], Dict[str, List[Psm]]): The peptide spectrum matches
        num_threads (int, optional): The number of threads to use. Defaults to 4.

    Returns:
        List[str]: The list of peptide sequences
    """

    psms = []

    if isinstance(psm_collection, dict):
        for _, psm_candidates in psm_collection.items():
            psms.extend(psm_candidates)

    else:
        psms = psm_collection

    return psc.get_psm_sequences_par([psm.get_py_ptr() for psm in psms], num_threads)

def get_spec_idx(psm_collection: Union[List[Psm], Dict[str, List[Psm]]], num_threads: int = 4) -> List[str]:
    """Get the spectrum indices from a list of peptide spectrum matches

    Args:
        psm_collection (Union[List[Psm], Dict[str, List[Psm]]): The peptide spectrum matches
        num_threads (int, optional): The number of threads to use. Defaults to 4.

    Returns:
        List[str]: The list of spectrum indices
    """

    psms = []

    if isinstance(psm_collection, dict):
        for _, psm_candidates in psm_collection.items():
            psms.extend(psm_candidates)

    else:
        psms = psm_collection

    return psc.get_psm_spec_idx_par([psm.get_py_ptr() for psm in psms], num_threads)


def psm_collection_to_pandas(psm_collection: Union[List[Psm], Dict[str, List[Psm]]],
                             num_threads: int = 4) -> pd.DataFrame:
    """Convert a list of peptide spectrum matches to a pandas dataframe

    Args:
        psm_collection (Union[List[Psm], Dict[str, List[Psm]]): The peptide spectrum matches
        num_threads (int, optional): The number of threads to use. Defaults to 4.

    Returns:
        pd.DataFrame: The pandas dataframe
    """

    psms = []

    if isinstance(psm_collection, dict):
        for _, item in psm_collection.items():
            psms.extend(item)

    else:
        psms = psm_collection

    # extract the numeric features
    D = psc.psms_to_feature_matrix([psm.get_py_ptr() for psm in psms], num_threads=num_threads)

    # extract the peptide sequences and spectrum indices
    sequence = psc.get_psm_sequences_par([psm.get_py_ptr() for psm in psms], num_threads=num_threads)
    sequence_modified = psc.get_psm_sequences_modified_par([psm.get_py_ptr() for psm in psms], num_threads=num_threads)
    sequence_decoy = psc.get_psm_sequences_decoy_par([psm.get_py_ptr() for psm in psms], num_threads=num_threads)
    sequence_decoy_modified = psc.get_psm_sequences_decoy_modified_par([psm.get_py_ptr() for psm in psms], num_threads=num_threads)
    spec_idx = psc.get_psm_spec_idx_par([psm.get_py_ptr() for psm in psms], num_threads=num_threads)
    proteins = psc.get_psm_proteins_par([psm.get_py_ptr() for psm in psms], num_threads=num_threads)

    # get the feature names
    names = psms[0].get_feature_names()

    # create the pandas dataframe
    PSM_pandas = pd.DataFrame(D, columns=names)

    # convert the decoy column to boolean
    PSM_pandas["decoy"] = [True if d == -1 else False for d in PSM_pandas.decoy]

    # add the sequence and spectrum index columns
    PSM_pandas.insert(0, "spec_idx", spec_idx)
    PSM_pandas.insert(1, "match_idx", sequence)
    PSM_pandas.insert(2,"match_identity_candidates", proteins)
    PSM_pandas.insert(3, "sequence", sequence)
    PSM_pandas.insert(4, "sequence_modified", sequence_modified)
    PSM_pandas.insert(5, "sequence_decoy", sequence_decoy)
    PSM_pandas.insert(6, "sequence_decoy_modified", sequence_decoy_modified)
    PSM_pandas.insert(7, "proteins", proteins)

    return PSM_pandas

def split_fasta(fasta: str, num_splits: int = 16, randomize: bool = True, verbose: bool = False) -> List[str]:
    """ Split a fasta file into multiple fasta files.
    Args:
        fasta: Fasta file as string.
        num_splits: Number of splits fasta file should be split into.
        randomize: Whether to randomize the order of sequences before splitting.

    Returns:
        List of fasta files as strings, will contain num_splits fasta files with equal number of sequences.
    """

    if num_splits == 1:
        return [fasta]

    split_strings = re.split(r'\n>', fasta)

    if verbose:
        print(f"Total number of sequences: {len(split_strings)} ...")

    if randomize:
        np.random.shuffle(split_strings)

    if not split_strings[0].startswith('>'):
        split_strings[0] = '>' + split_strings[0]

    total_items = len(split_strings)
    items_per_batch = total_items // num_splits
    remainder = total_items % num_splits

    fastas = []
    start_index = 0

    for i in range(num_splits):
        extra = 1 if i < remainder else 0
        stop_index = start_index + items_per_batch + extra

        if start_index >= total_items:
            break

        batch = '\n>'.join(split_strings[start_index:stop_index])

        if not batch.startswith('>'):
            batch = '>' + batch

        fastas.append(batch)
        start_index = stop_index

    return fastas

def generate_search_configurations(
    fasta_path: str,
    num_splits: int = 25,
    missed_cleavages: int = 1,
    min_len: int = 8,
    max_len: int = 30,
    cleave_at: str = "KR",
    restrict: str = "P",
    c_terminal: bool = True,
    static_mods: dict = {"C": "[UNIMOD:4]"},  # Static cysteine modification
    variable_mods: dict = {"M": ["[UNIMOD:35]"]},  # Oxidation on methionine
    bucket_size: int = 2**14,
    generate_decoys: bool = True,
    randomize_split: bool = True,
) -> Iterator[SageSearchConfiguration]:
    """
    Generates an iterator of indexed databases for each split of a FASTA file.

    Args:
        fasta_path (str): Path to the unsplit FASTA file.
        num_splits (int): Number of splits for the FASTA file. Default is 25.
        missed_cleavages (int): Allowed missed cleavages for the enzyme. Default is 2.
        min_len (int): Minimum peptide length. Default is 5.
        max_len (int): Maximum peptide length. Default is 50.
        cleave_at (str): Amino acids where cleavage occurs. Default is "KR".
        restrict (str): Restriction amino acids. Default is "P".
        c_terminal (bool): Whether cleavage is C-terminal. Default is True.
        static_mods (dict): Static modifications in UNIMOD notation. Default is {"C": "[UNIMOD:4]"}.
        variable_mods (dict): Variable modifications in UNIMOD notation. Default is {"M": ["[UNIMOD:35]"]}.
        bucket_size (int): Size of the bucket for indexing. Default is 2^14.
        generate_decoys (bool): Whether to generate decoys in the database. Default is True.
        randomize_split (bool): Whether to randomize the order of sequences before splitting. Default is True.

    Yields:
        SageSearchConfiguration: Indexed database configuration for each split.
    """
    # Read the full FASTA file
    with open(fasta_path, "r") as infile:
        fasta = infile.read()

    # Split the FASTA file
    fastas = split_fasta(fasta, num_splits=num_splits, randomize=randomize_split)

    # Configure the enzyme builder
    enzyme_builder = EnzymeBuilder(
        missed_cleavages=missed_cleavages,
        min_len=min_len,
        max_len=max_len,
        cleave_at=cleave_at,
        restrict=restrict,
        c_terminal=c_terminal,
    )

    # Generate configurations for each split
    for fasta in fastas:
        sage_config = SageSearchConfiguration(
            fasta=fasta,
            static_mods=static_mods,
            variable_mods=variable_mods,
            enzyme_builder=enzyme_builder,
            generate_decoys=generate_decoys,
            bucket_size=bucket_size,
        )

        # Generate the indexed database
        indexed_db = sage_config.generate_indexed_database()

        # Yield the configuration with the indexed database
        yield indexed_db


def compress_psms(psms: List[Psm]) -> bytes:
    """Compress a list of peptide spectrum matches.

    Args:
        psms (List[Psm]): The peptide spectrum matches

    Returns:
        List[Psm]: The compressed peptide spectrum matches
    """
    return psc.py_compress_psms([psm.get_py_ptr() for psm in psms])

def decompress_psms(data: bytes) -> List[Psm]:
    """Decompress a list of peptide spectrum matches.

    Args:
        data (bytes): The compressed peptide spectrum matches

    Returns:
        List[Psm]: The decompressed peptide spectrum matches
    """
    return [Psm.from_py_ptr(p) for p in psc.py_decompress_psms(data)]
