"""
Fragment ion intensity prediction interface layer.

This module provides the interface between sagepy and external intensity prediction
libraries (like imspy). It handles:
- Exporting peptide sequences for prediction (PredictionRequest)
- Importing predicted intensities (PredictionResult)
- Serializing predictions to the .sagi binary format
"""

from __future__ import annotations

import struct
from collections import defaultdict
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from .database import IndexedDatabase


# Ion kind codes matching SAGE convention
ION_KIND_A = 0
ION_KIND_B = 1
ION_KIND_C = 2
ION_KIND_X = 3
ION_KIND_Y = 4
ION_KIND_Z = 5

DEFAULT_ION_KINDS = [ION_KIND_B, ION_KIND_Y]


@dataclass
class PredictionRequest:
    """
    Batch of peptides for intensity prediction.

    Attributes
    ----------
    sequences : np.ndarray
        Modified peptide sequences with UNIMOD notation (e.g., "PEPTC[UNIMOD:4]IDEK")
    charges : np.ndarray
        Precursor charge states (parallel to sequences)
    peptide_indices : np.ndarray
        IndexedDatabase peptide indices for mapping predictions back
    """

    sequences: np.ndarray
    charges: np.ndarray
    peptide_indices: np.ndarray

    def __post_init__(self):
        self.sequences = np.asarray(self.sequences)
        self.charges = np.asarray(self.charges, dtype=np.int32)
        self.peptide_indices = np.asarray(self.peptide_indices, dtype=np.int64)

        if not (len(self.sequences) == len(self.charges) == len(self.peptide_indices)):
            raise ValueError(
                f"Array length mismatch: sequences={len(self.sequences)}, "
                f"charges={len(self.charges)}, peptide_indices={len(self.peptide_indices)}"
            )

    def to_dataframe(self) -> pd.DataFrame:
        """Convert to DataFrame for easy inspection/manipulation."""
        return pd.DataFrame(
            {
                "sequence": self.sequences,
                "charge": self.charges,
                "peptide_idx": self.peptide_indices,
            }
        )

    @classmethod
    def from_dataframe(cls, df: pd.DataFrame) -> PredictionRequest:
        """Create PredictionRequest from DataFrame."""
        return cls(
            sequences=df["sequence"].values,
            charges=df["charge"].values,
            peptide_indices=df["peptide_idx"].values,
        )

    def __len__(self) -> int:
        return len(self.sequences)

    def __repr__(self) -> str:
        return (
            f"PredictionRequest(n_entries={len(self)}, "
            f"unique_peptides={len(np.unique(self.peptide_indices))}, "
            f"charges={sorted(np.unique(self.charges).tolist())})"
        )


@dataclass
class PredictionResult:
    """
    Predicted intensities for a batch of peptides.

    Attributes
    ----------
    peptide_indices : np.ndarray
        IndexedDatabase peptide indices (must match request)
    charges : np.ndarray
        Precursor charges (must match request)
    intensities : list[np.ndarray]
        Fragment intensities, one array per peptide.
        Each array has shape [n_ion_kinds, seq_len-1, max_frag_charge]
        Values should be normalized (0.0 to 1.0)
    ion_kinds : list[int]
        Ion type codes used (default: [1, 4] for B and Y ions)
    max_fragment_charge : int
        Maximum fragment charge state predicted
    """

    peptide_indices: np.ndarray
    charges: np.ndarray
    intensities: list[np.ndarray]
    ion_kinds: list[int] = field(default_factory=lambda: DEFAULT_ION_KINDS.copy())
    max_fragment_charge: int = 2

    def __post_init__(self):
        self.peptide_indices = np.asarray(self.peptide_indices, dtype=np.int64)
        self.charges = np.asarray(self.charges, dtype=np.int32)

        if not (len(self.peptide_indices) == len(self.charges) == len(self.intensities)):
            raise ValueError(
                f"Array length mismatch: peptide_indices={len(self.peptide_indices)}, "
                f"charges={len(self.charges)}, intensities={len(self.intensities)}"
            )

    def __len__(self) -> int:
        return len(self.intensities)

    def __repr__(self) -> str:
        return (
            f"PredictionResult(n_entries={len(self)}, "
            f"ion_kinds={self.ion_kinds}, "
            f"max_fragment_charge={self.max_fragment_charge})"
        )


@dataclass
class PredictionResultV2:
    """
    V2 Predicted intensities keyed by (unimod_sequence, charge) instead of peptide_idx.

    This format decouples predictions from the database order, allowing:
    - Database chunking (splitting peptides into batches)
    - Reusing predictions across different database configurations
    - Flexible parallel processing workflows
    - **Different intensity predictions for different PTM states**

    Attributes
    ----------
    predictions : dict[tuple[str, int], np.ndarray]
        Map from (unimod_sequence, precursor_charge) to intensity tensor.
        Each tensor has shape [n_ion_kinds, seq_len-1, max_frag_charge].
        The unimod_sequence should include UNIMOD notation for modifications,
        e.g., "PEPTC[UNIMOD:4]IDEK" for carbamidomethylated cysteine.
        This allows different modifications to have different intensity predictions.
    ion_kinds : list[int]
        Ion type codes used (default: [1, 4] for B and Y ions)
    max_fragment_charge : int
        Maximum fragment charge state predicted

    Example
    -------
    >>> predictions = {
    ...     ("PEPTIDEK", 2): intensity_unmodified,
    ...     ("PEPTC[UNIMOD:4]IDEK", 2): intensity_carbamidomethyl,
    ...     ("PEPTM[UNIMOD:35]IDEK", 2): intensity_oxidized,
    ... }
    >>> result = PredictionResultV2(predictions=predictions)
    """

    predictions: dict[tuple[str, int], np.ndarray]
    ion_kinds: list[int] = field(default_factory=lambda: DEFAULT_ION_KINDS.copy())
    max_fragment_charge: int = 2

    def __len__(self) -> int:
        return len(self.predictions)

    def __repr__(self) -> str:
        return (
            f"PredictionResultV2(n_entries={len(self)}, "
            f"ion_kinds={self.ion_kinds}, "
            f"max_fragment_charge={self.max_fragment_charge})"
        )


def compute_key_hash(sequence: str, charge: int) -> int:
    """
    Compute hash for (sequence, charge) key matching Rust implementation.

    This uses the same hashing algorithm as the Rust implementation,
    ensuring consistency between Python and Rust lookups.

    Parameters
    ----------
    sequence : str
        Raw peptide sequence (NOT UNIMOD notation)
    charge : int
        Precursor charge state

    Returns
    -------
    int
        Hash value for the (sequence, charge) key
    """
    import sagepy_connector as psc

    return psc.py_intensity.PyPredictedIntensityStore.compute_key_hash(sequence, charge)


def create_prediction_request(
    indexed_db: IndexedDatabase,
    charges: list[int] | None = None,
    include_decoys: bool = False,
) -> PredictionRequest:
    """
    Create prediction request from IndexedDatabase.

    Parameters
    ----------
    indexed_db : IndexedDatabase
        The peptide database
    charges : list[int], optional
        Precursor charges to predict for each peptide. Default: [2, 3]
    include_decoys : bool
        Whether to include decoy sequences. Default: False

    Returns
    -------
    PredictionRequest
        Batch of (sequence, charge, peptide_idx) for prediction
    """
    if charges is None:
        charges = [2, 3]

    sequences = indexed_db.peptides_as_string()
    n_peptides = len(sequences)

    if not include_decoys:
        decoy_tag = indexed_db.decoy_tag
        # Create mask for non-decoy sequences
        mask = np.array([decoy_tag not in seq for seq in sequences])
        sequences = sequences[mask]
        indices = np.arange(n_peptides)[mask]
    else:
        indices = np.arange(n_peptides)

    # Expand for each charge state
    n_charges = len(charges)
    expanded_sequences = np.repeat(sequences, n_charges)
    expanded_charges = np.tile(charges, len(sequences))
    expanded_indices = np.repeat(indices, n_charges)

    return PredictionRequest(
        sequences=expanded_sequences,
        charges=expanded_charges,
        peptide_indices=expanded_indices,
    )


def get_sequence_length(sequence: str) -> int:
    """
    Get the number of amino acids in a modified sequence.

    Handles UNIMOD notation like "PEPTC[UNIMOD:4]IDEK" -> 9 AAs.
    """
    length = 0
    i = 0
    while i < len(sequence):
        if sequence[i] == "[":
            # Skip modification block
            while i < len(sequence) and sequence[i] != "]":
                i += 1
            i += 1  # Skip closing bracket
        elif sequence[i].isupper():
            length += 1
            i += 1
        else:
            i += 1
    return length


def validate_prediction_result(
    request: PredictionRequest, result: PredictionResult
) -> bool:
    """
    Validate that result matches request.

    Parameters
    ----------
    request : PredictionRequest
        The original prediction request
    result : PredictionResult
        The prediction result to validate

    Returns
    -------
    bool
        True if validation passes

    Raises
    ------
    ValueError
        If validation fails with description of the mismatch
    """
    if len(request) != len(result.intensities):
        raise ValueError(
            f"Length mismatch: request has {len(request)}, "
            f"result has {len(result.intensities)}"
        )

    if not np.array_equal(request.peptide_indices, result.peptide_indices):
        raise ValueError("Peptide indices don't match")

    if not np.array_equal(request.charges, result.charges):
        raise ValueError("Charges don't match")

    # Validate intensity shapes
    n_ion_kinds = len(result.ion_kinds)
    max_frag_charge = result.max_fragment_charge

    for i, (seq, intensity) in enumerate(zip(request.sequences, result.intensities)):
        seq_len = get_sequence_length(seq)
        expected_shape = (n_ion_kinds, seq_len - 1, max_frag_charge)

        if intensity.shape != expected_shape:
            raise ValueError(
                f"Intensity shape mismatch at index {i} (sequence '{seq}'): "
                f"expected {expected_shape}, got {intensity.shape}"
            )

    return True


def aggregate_predictions_by_peptide(
    result: PredictionResult,
    aggregation: Literal["max_charge", "min_charge", "mean"] = "max_charge",
) -> dict[int, np.ndarray]:
    """
    Aggregate predictions to one per peptide_idx.

    Since the .sagi format stores one prediction per peptide (not per charge),
    this function aggregates multiple charge predictions for the same peptide.

    Parameters
    ----------
    result : PredictionResult
        Prediction result with potentially multiple entries per peptide
    aggregation : str
        Aggregation method:
        - 'max_charge': keep prediction for highest precursor charge
        - 'min_charge': keep prediction for lowest precursor charge
        - 'mean': average across all charge predictions

    Returns
    -------
    dict[int, np.ndarray]
        Mapping from peptide_idx to aggregated intensity array
    """
    grouped: dict[int, list[tuple[int, np.ndarray]]] = defaultdict(list)

    for idx, charge, intensity in zip(
        result.peptide_indices, result.charges, result.intensities
    ):
        grouped[idx].append((charge, intensity))

    aggregated: dict[int, np.ndarray] = {}

    for idx, predictions in grouped.items():
        if aggregation == "max_charge":
            max_charge = max(p[0] for p in predictions)
            aggregated[idx] = next(p[1] for p in predictions if p[0] == max_charge)
        elif aggregation == "min_charge":
            min_charge = min(p[0] for p in predictions)
            aggregated[idx] = next(p[1] for p in predictions if p[0] == min_charge)
        elif aggregation == "mean":
            aggregated[idx] = np.mean([p[1] for p in predictions], axis=0)
        else:
            raise ValueError(f"Unknown aggregation method: {aggregation}")

    return aggregated


def write_intensity_file(
    path: str,
    predictions: list[np.ndarray],
    peptide_lengths: list[int],
    max_charge: int = 2,
    ion_kinds: list[int] | None = None,
) -> None:
    """
    Write predicted intensities to Sage binary format (.sagi).

    Parameters
    ----------
    path : str
        Output file path (typically .sagi extension)
    predictions : list[np.ndarray]
        One array per peptide, shape [n_ion_kinds, peptide_len-1, max_charge].
        Values should be normalized intensities (e.g., 0.0 to 1.0).
        Must be in peptide index order (matching IndexedDatabase).
    peptide_lengths : list[int]
        Length of each peptide (for validation)
    max_charge : int
        Maximum fragment charge state. Default: 2
    ion_kinds : list[int], optional
        Ion type codes. Default: [1, 4] for B and Y ions

    Example
    -------
    >>> # For 3 peptides with B and Y ions, max charge 2
    >>> predictions = [
    ...     np.random.rand(2, len1-1, 2).astype(np.float32),
    ...     np.random.rand(2, len2-1, 2).astype(np.float32),
    ...     np.random.rand(2, len3-1, 2).astype(np.float32),
    ... ]
    >>> write_intensity_file("predictions.sagi", predictions, [len1, len2, len3])
    """
    if ion_kinds is None:
        ion_kinds = DEFAULT_ION_KINDS.copy()

    if len(predictions) != len(peptide_lengths):
        raise ValueError(
            f"Length mismatch: {len(predictions)} predictions, "
            f"{len(peptide_lengths)} peptide lengths"
        )

    with open(path, "wb") as f:
        # Header
        f.write(struct.pack("<I", 0x49474153))  # magic "SAGI"
        f.write(struct.pack("<I", 1))  # version
        f.write(struct.pack("<Q", len(predictions)))  # peptide_count
        f.write(struct.pack("<B", max_charge))
        f.write(struct.pack("<B", len(ion_kinds)))
        for k in ion_kinds:
            f.write(struct.pack("<B", k))

        # Calculate and write offsets
        offsets = []
        current_offset = 0
        for pred in predictions:
            offsets.append(current_offset)
            current_offset += pred.size * 4  # f32 = 4 bytes

        for off in offsets:
            f.write(struct.pack("<Q", off))

        # Write data
        for i, pred in enumerate(predictions):
            # Ensure correct memory layout and dtype
            pred_f32 = pred.astype("<f4")
            f.write(pred_f32.tobytes())


def read_intensity_file(path: str) -> dict:
    """
    Read predicted intensities from Sage binary format (.sagi).

    Auto-detects the file format version (V1 or V2) and reads accordingly.

    Parameters
    ----------
    path : str
        Path to .sagi file

    Returns
    -------
    dict
        For V1 files:
        - version: int (1)
        - peptide_count: int
        - max_charge: int
        - ion_kinds: list[int]
        - offsets: np.ndarray of uint64
        - data: np.ndarray of float32
        - is_key_based: bool (False)

        For V2 files:
        - version: int (2)
        - entry_count: int
        - max_charge: int
        - ion_kinds: list[int]
        - is_key_based: bool (True)

    Notes
    -----
    For V2 files, actual intensity data is accessed via the PredictedIntensityStore
    class using key-based lookups (get_intensity_by_key). This function only reads
    the header metadata.
    """
    with open(path, "rb") as f:
        magic = struct.unpack("<I", f.read(4))[0]
        if magic != 0x49474153:
            raise ValueError(f"Invalid magic: {hex(magic)}, expected 0x49474153 ('SAGI')")

        version = struct.unpack("<I", f.read(4))[0]

        if version == 1:
            # V1 format: positional indexing
            peptide_count = struct.unpack("<Q", f.read(8))[0]
            max_charge = struct.unpack("<B", f.read(1))[0]
            ion_kind_count = struct.unpack("<B", f.read(1))[0]
            ion_kinds = [struct.unpack("<B", f.read(1))[0] for _ in range(ion_kind_count)]

            offsets = np.frombuffer(f.read(peptide_count * 8), dtype="<u8")
            data = np.frombuffer(f.read(), dtype="<f4")

            return {
                "version": 1,
                "peptide_count": peptide_count,
                "max_charge": max_charge,
                "ion_kinds": ion_kinds,
                "offsets": offsets,
                "data": data,
                "is_key_based": False,
            }

        elif version == 2:
            # V2 format: key-based indexing
            entry_count = struct.unpack("<Q", f.read(8))[0]
            max_charge = struct.unpack("<B", f.read(1))[0]
            ion_kind_count = struct.unpack("<B", f.read(1))[0]
            ion_kinds = [struct.unpack("<B", f.read(1))[0] for _ in range(ion_kind_count)]

            return {
                "version": 2,
                "entry_count": entry_count,
                "max_charge": max_charge,
                "ion_kinds": ion_kinds,
                "is_key_based": True,
            }

        else:
            raise ValueError(f"Unsupported version: {version}")


def write_intensity_file_v2(
    path: str,
    predictions: dict[tuple[str, int], np.ndarray],
    max_charge: int = 2,
    ion_kinds: list[int] | None = None,
) -> None:
    """
    Write predicted intensities to Sage V2 binary format (.sagi).

    V2 format uses (unimod_sequence, charge) as the key, allowing database chunking,
    reuse of predictions across different database configurations, and importantly,
    **different intensity predictions for different PTM states**.

    Parameters
    ----------
    path : str
        Output file path (typically .sagi extension)
    predictions : dict[tuple[str, int], np.ndarray]
        Map from (unimod_sequence, precursor_charge) to intensity tensor.
        Each tensor has shape [n_ion_kinds, seq_len-1, max_charge].
        The unimod_sequence should include UNIMOD notation for modifications,
        e.g., "PEPTC[UNIMOD:4]IDEK" for carbamidomethylated cysteine.
    max_charge : int
        Maximum fragment charge state. Default: 2
    ion_kinds : list[int], optional
        Ion type codes. Default: [1, 4] for B and Y ions

    Example
    -------
    >>> predictions = {
    ...     ("PEPTIDEK", 2): intensity_unmodified,
    ...     ("PEPTC[UNIMOD:4]IDEK", 2): intensity_carbamidomethyl,
    ...     ("PEPTM[UNIMOD:35]IDEK", 2): intensity_oxidized,
    ... }
    >>> write_intensity_file_v2("predictions_v2.sagi", predictions)
    """
    if ion_kinds is None:
        ion_kinds = DEFAULT_ION_KINDS.copy()

    # Convert predictions to format expected by Rust
    # entries: list of (sequence_bytes, charge, peptide_len, flat_intensities)
    entries = []
    for (sequence, charge), intensity in predictions.items():
        sequence_bytes = sequence.encode("utf-8") if isinstance(sequence, str) else sequence
        peptide_len = intensity.shape[1] + 1  # intensity has shape [ion_kinds, seq_len-1, max_charge]
        flat_intensities = intensity.astype(np.float32).flatten().tolist()
        entries.append((list(sequence_bytes), charge, peptide_len, flat_intensities))

    # Use Rust implementation to create and write the store
    import sagepy_connector as psc

    store = psc.py_intensity.PyPredictedIntensityStore.from_raw_predictions_v2(
        entries, max_charge, ion_kinds
    )
    store.write(path)


def read_intensity_file_v2(path: str) -> dict:
    """
    Read predicted intensities from Sage V2 binary format (.sagi).

    Parameters
    ----------
    path : str
        Path to .sagi V2 file

    Returns
    -------
    dict
        Dictionary with keys:
        - version: int (2 for V2 format)
        - entry_count: int
        - max_charge: int
        - ion_kinds: list[int]
        - is_key_based: bool (True for V2)

    Notes
    -----
    For V2 files, actual intensity data is accessed via the PredictedIntensityStore
    class using key-based lookups (get_intensity_by_key).
    """
    with open(path, "rb") as f:
        magic = struct.unpack("<I", f.read(4))[0]
        if magic != 0x49474153:
            raise ValueError(f"Invalid magic: {hex(magic)}, expected 0x49474153 ('SAGI')")

        version = struct.unpack("<I", f.read(4))[0]
        if version != 2:
            raise ValueError(f"Expected V2 format (version 2), got version {version}")

        entry_count = struct.unpack("<Q", f.read(8))[0]
        max_charge = struct.unpack("<B", f.read(1))[0]
        ion_kind_count = struct.unpack("<B", f.read(1))[0]
        ion_kinds = [struct.unpack("<B", f.read(1))[0] for _ in range(ion_kind_count)]

        return {
            "version": version,
            "entry_count": entry_count,
            "max_charge": max_charge,
            "ion_kinds": ion_kinds,
            "is_key_based": True,
        }


def get_peptide_intensities(
    intensity_data: dict,
    peptide_idx: int,
    peptide_length: int,
) -> np.ndarray:
    """
    Extract intensity array for a specific peptide from loaded .sagi data.

    Parameters
    ----------
    intensity_data : dict
        Data returned by read_intensity_file()
    peptide_idx : int
        Index of the peptide
    peptide_length : int
        Length of the peptide sequence

    Returns
    -------
    np.ndarray
        Intensity array with shape [n_ion_kinds, peptide_len-1, max_charge]
    """
    n_ion_kinds = len(intensity_data["ion_kinds"])
    max_charge = intensity_data["max_charge"]
    n_positions = peptide_length - 1

    offset = int(intensity_data["offsets"][peptide_idx])
    n_values = n_ion_kinds * n_positions * max_charge

    # Offset is in bytes, data is float32 (4 bytes)
    start_idx = offset // 4
    end_idx = start_idx + n_values

    flat_data = intensity_data["data"][start_idx:end_idx]
    return flat_data.reshape(n_ion_kinds, n_positions, max_charge)


def create_uniform_intensity(
    peptide_length: int,
    ion_kinds: list[int] | None = None,
    max_charge: int = 2,
    value: float = 1.0,
) -> np.ndarray:
    """
    Create uniform intensity array for a peptide.

    Useful as fallback for peptides without predictions.

    Parameters
    ----------
    peptide_length : int
        Length of the peptide sequence
    ion_kinds : list[int], optional
        Ion type codes. Default: [1, 4] for B and Y
    max_charge : int
        Maximum fragment charge. Default: 2
    value : float
        Uniform intensity value. Default: 1.0

    Returns
    -------
    np.ndarray
        Uniform intensity array with shape [n_ion_kinds, peptide_len-1, max_charge]
    """
    if ion_kinds is None:
        ion_kinds = DEFAULT_ION_KINDS.copy()

    n_ion_kinds = len(ion_kinds)
    n_positions = peptide_length - 1

    return np.full(
        (n_ion_kinds, n_positions, max_charge), value, dtype=np.float32
    )


def build_intensity_file_from_result(
    indexed_db: IndexedDatabase,
    result: PredictionResult,
    output_path: str,
    aggregation: Literal["max_charge", "min_charge", "mean"] = "max_charge",
) -> None:
    """
    Build a complete .sagi file from prediction results.

    This is a convenience function that:
    1. Aggregates predictions by peptide index
    2. Fills in uniform intensities for missing peptides (decoys)
    3. Writes the .sagi file

    Parameters
    ----------
    indexed_db : IndexedDatabase
        The peptide database (needed for peptide count and lengths)
    result : PredictionResult
        Prediction results from external predictor
    output_path : str
        Path for output .sagi file
    aggregation : str
        How to aggregate multiple charge predictions. Default: 'max_charge'
    """
    aggregated = aggregate_predictions_by_peptide(result, aggregation)

    sequences = indexed_db.peptides_as_string()
    n_peptides = len(sequences)

    predictions = []
    peptide_lengths = []

    for idx in range(n_peptides):
        seq_len = get_sequence_length(sequences[idx])
        peptide_lengths.append(seq_len)

        if idx in aggregated:
            predictions.append(aggregated[idx])
        else:
            # Missing prediction (likely decoy) - use uniform
            predictions.append(
                create_uniform_intensity(
                    seq_len,
                    ion_kinds=result.ion_kinds,
                    max_charge=result.max_fragment_charge,
                )
            )

    write_intensity_file(
        output_path,
        predictions,
        peptide_lengths,
        max_charge=result.max_fragment_charge,
        ion_kinds=result.ion_kinds,
    )


class PredictedIntensityStore:
    """
    Store for predicted fragment ion intensities.

    Wraps the Rust IntensityStore loaded from a .sagi file.
    Supports both V1 (positional indexing) and V2 (key-based indexing) formats.
    This store is used during scoring to weight matched fragment peaks
    by their predicted intensities.

    Example (V1 - positional lookup)
    --------------------------------
    >>> store = PredictedIntensityStore("predictions.sagi")
    >>> print(f"Loaded {store.peptide_count} peptides")
    >>> intensity = store.get_intensity(
    ...     peptide_idx=0,
    ...     peptide_len=10,
    ...     ion_kind=ION_KIND_B,
    ...     position=3,
    ...     charge=1
    ... )

    Example (V2 - key-based lookup)
    -------------------------------
    >>> store = PredictedIntensityStore("predictions_v2.sagi")
    >>> print(f"Is key-based: {store.is_key_based}")
    >>> intensity = store.get_intensity_by_key(
    ...     sequence="PEPTIDEK",
    ...     precursor_charge=2,
    ...     ion_kind=ION_KIND_B,
    ...     position=3,
    ...     fragment_charge=1
    ... )
    """

    def __init__(self, path: str):
        """
        Load intensity store from a .sagi binary file.

        Auto-detects the file format version (V1 or V2) and loads accordingly.

        Parameters
        ----------
        path : str
            Path to the .sagi file created by write_intensity_file(),
            write_intensity_file_v2(), or from_predictions_v2()
        """
        import sagepy_connector as psc

        self.__py_ptr = psc.py_intensity.PyPredictedIntensityStore.load(path)

    @classmethod
    def from_py_ptr(cls, py_ptr) -> "PredictedIntensityStore":
        """Create a PredictedIntensityStore from an existing PyO3 pointer."""
        instance = cls.__new__(cls)
        instance.__py_ptr = py_ptr
        return instance

    @classmethod
    def uniform(
        cls,
        peptide_lengths: list[int],
        max_charge: int = 3,
        ion_kinds: list[int] | None = None,
    ) -> "PredictedIntensityStore":
        """
        Create a uniform V1 intensity store where all intensities are 1.0.

        This is useful for testing the weighted scoring code path without
        actual predictions. With uniform intensities, weighted scores should
        be identical to unweighted scores.

        Parameters
        ----------
        peptide_lengths : list[int]
            Length of each peptide in the database. Can be obtained from
            IndexedDatabase via `[len(seq) for seq in db.peptides_as_string()]`
        max_charge : int
            Maximum fragment charge state (default: 3)
        ion_kinds : list[int], optional
            Ion type codes (default: [1, 4] for B and Y ions)

        Returns
        -------
        PredictedIntensityStore
            V1 store with all intensities set to 1.0

        Example
        -------
        >>> # Create uniform store for testing
        >>> peptide_lengths = [len(seq) for seq in indexed_db.peptides_as_string()]
        >>> store = PredictedIntensityStore.uniform(peptide_lengths)
        >>> # Use with weighted scoring - results should match unweighted
        >>> scorer = Scorer(score_type=ScoreType("weightedhyperscore"))
        >>> features = scorer.score(db, spectrum, intensity_store=store)
        """
        import sagepy_connector as psc

        py_ptr = psc.py_intensity.PyPredictedIntensityStore.uniform(
            peptide_lengths, max_charge, ion_kinds
        )
        return cls.from_py_ptr(py_ptr)

    @classmethod
    def from_predictions_v2(
        cls,
        predictions: dict[tuple[str, int], np.ndarray],
        max_charge: int = 2,
        ion_kinds: list[int] | None = None,
    ) -> "PredictedIntensityStore":
        """
        Create a V2 intensity store from prediction dict.

        V2 format uses (unimod_sequence, charge) as the key, allowing database chunking,
        reuse of predictions across different database configurations, and importantly,
        **different intensity predictions for different PTM states**.

        Parameters
        ----------
        predictions : dict[tuple[str, int], np.ndarray]
            Map from (unimod_sequence, precursor_charge) to intensity tensor.
            Each tensor has shape [n_ion_kinds, seq_len-1, max_charge].
            The unimod_sequence should include UNIMOD notation for modifications,
            e.g., "PEPTC[UNIMOD:4]IDEK" for carbamidomethylated cysteine.
        max_charge : int
            Maximum fragment charge state (default: 2)
        ion_kinds : list[int], optional
            Ion type codes (default: [1, 4] for B and Y ions)

        Returns
        -------
        PredictedIntensityStore
            V2 store with key-based indexing

        Example
        -------
        >>> # Different predictions for different modification states
        >>> predictions = {
        ...     ("PEPTIDEK", 2): intensity_unmodified,
        ...     ("PEPTC[UNIMOD:4]IDEK", 2): intensity_carbamidomethyl,
        ...     ("PEPTM[UNIMOD:35]IDEK", 2): intensity_oxidized,
        ... }
        >>> store = PredictedIntensityStore.from_predictions_v2(predictions)
        >>> # Use with different database chunks
        >>> for chunk in database_chunks:
        ...     scores = scorer.score(chunk, spectra, intensity_store=store)
        """
        import sagepy_connector as psc

        if ion_kinds is None:
            ion_kinds = DEFAULT_ION_KINDS.copy()

        # Convert predictions to format expected by Rust
        entries = []
        for (sequence, charge), intensity in predictions.items():
            sequence_bytes = (
                sequence.encode("utf-8") if isinstance(sequence, str) else sequence
            )
            peptide_len = intensity.shape[1] + 1
            flat_intensities = intensity.astype(np.float32).flatten().tolist()
            entries.append((list(sequence_bytes), charge, peptide_len, flat_intensities))

        py_ptr = psc.py_intensity.PyPredictedIntensityStore.from_raw_predictions_v2(
            entries, max_charge, ion_kinds
        )
        return cls.from_py_ptr(py_ptr)

    def get_py_ptr(self):
        """Get the underlying PyO3 pointer for passing to Rust functions."""
        return self.__py_ptr

    @property
    def is_key_based(self) -> bool:
        """True if this is a V2 store with key-based indexing."""
        return self.__py_ptr.is_key_based

    def get_intensity(
        self,
        peptide_idx: int,
        peptide_len: int,
        ion_kind: int,
        position: int,
        charge: int,
    ) -> float | None:
        """
        Get predicted intensity for a specific fragment (V1 positional lookup).

        Parameters
        ----------
        peptide_idx : int
            Index of the peptide in IndexedDatabase
        peptide_len : int
            Length of the peptide sequence
        ion_kind : int
            Ion type code (0=A, 1=B, 2=C, 3=X, 4=Y, 5=Z)
        position : int
            Fragment position (0-indexed)
        charge : int
            Charge state (1-indexed)

        Returns
        -------
        float or None
            The predicted intensity, or None if not found or if this is a V2 store
        """
        return self.__py_ptr.get_intensity(
            peptide_idx, peptide_len, ion_kind, position, charge
        )

    def get_intensity_or_default(
        self,
        peptide_idx: int,
        peptide_len: int,
        ion_kind: int,
        position: int,
        charge: int,
        default: float = 1.0,
    ) -> float:
        """
        Get predicted intensity with fallback to default (V1 positional lookup).

        Parameters
        ----------
        peptide_idx : int
            Index of the peptide in IndexedDatabase
        peptide_len : int
            Length of the peptide sequence
        ion_kind : int
            Ion type code (0=A, 1=B, 2=C, 3=X, 4=Y, 5=Z)
        position : int
            Fragment position (0-indexed)
        charge : int
            Charge state (1-indexed)
        default : float
            Value to return if intensity not found. Default: 1.0

        Returns
        -------
        float
            The predicted intensity, or default if not found
        """
        result = self.__py_ptr.get_intensity_or_default(
            peptide_idx, peptide_len, ion_kind, position, charge
        )
        return result if result is not None else default

    def get_intensity_by_key(
        self,
        sequence: str,
        precursor_charge: int,
        ion_kind: int,
        position: int,
        fragment_charge: int,
    ) -> float | None:
        """
        Get predicted intensity by sequence and charge (V2 key-based lookup).

        Parameters
        ----------
        sequence : str
            Raw peptide sequence (NOT UNIMOD notation)
        precursor_charge : int
            Precursor charge state
        ion_kind : int
            Ion type code (0=A, 1=B, 2=C, 3=X, 4=Y, 5=Z)
        position : int
            Fragment position (0-indexed)
        fragment_charge : int
            Fragment charge state (1-indexed)

        Returns
        -------
        float or None
            The predicted intensity, or None if not found or if this is a V1 store
        """
        return self.__py_ptr.get_intensity_by_key(
            sequence, precursor_charge, ion_kind, position, fragment_charge
        )

    def get_intensity_by_key_or_default(
        self,
        sequence: str,
        precursor_charge: int,
        ion_kind: int,
        position: int,
        fragment_charge: int,
        default: float = 1.0,
    ) -> float:
        """
        Get predicted intensity by key with fallback to default (V2 lookup).

        Parameters
        ----------
        sequence : str
            Raw peptide sequence (NOT UNIMOD notation)
        precursor_charge : int
            Precursor charge state
        ion_kind : int
            Ion type code (0=A, 1=B, 2=C, 3=X, 4=Y, 5=Z)
        position : int
            Fragment position (0-indexed)
        fragment_charge : int
            Fragment charge state (1-indexed)
        default : float
            Value to return if intensity not found. Default: 1.0

        Returns
        -------
        float
            The predicted intensity, or default if not found
        """
        result = self.__py_ptr.get_intensity_by_key_or_default(
            sequence, precursor_charge, ion_kind, position, fragment_charge
        )
        return result if result is not None else default

    def contains_key(self, sequence: str, precursor_charge: int) -> bool:
        """
        Check if a key exists in the V2 store.

        Parameters
        ----------
        sequence : str
            Raw peptide sequence (NOT UNIMOD notation)
        precursor_charge : int
            Precursor charge state

        Returns
        -------
        bool
            True if the key exists, False for V1 stores or missing keys
        """
        return self.__py_ptr.contains_key(sequence, precursor_charge)

    def write(self, path: str) -> None:
        """
        Write the store to a .sagi binary file.

        The file format version (V1 or V2) is determined by the store type.

        Parameters
        ----------
        path : str
            Output file path
        """
        self.__py_ptr.write(path)

    @property
    def entry_count(self) -> int:
        """Number of entries in the store (peptides for V1, keys for V2)."""
        return self.__py_ptr.entry_count

    @property
    def peptide_count(self) -> int:
        """Number of peptides in the store (V1 only, returns 0 for V2)."""
        return self.__py_ptr.peptide_count

    @property
    def max_charge(self) -> int:
        """Maximum fragment charge state stored."""
        return self.__py_ptr.max_charge

    @property
    def ion_kinds(self) -> list[int]:
        """Ion type codes stored (0=A, 1=B, 2=C, 3=X, 4=Y, 5=Z)."""
        return self.__py_ptr.ion_kinds

    def __repr__(self) -> str:
        if self.is_key_based:
            return (
                f"PredictedIntensityStore(v2, entry_count={self.entry_count}, "
                f"max_charge={self.max_charge}, ion_kinds={self.ion_kinds})"
            )
        else:
            return (
                f"PredictedIntensityStore(v1, peptide_count={self.peptide_count}, "
                f"max_charge={self.max_charge}, ion_kinds={self.ion_kinds})"
            )
