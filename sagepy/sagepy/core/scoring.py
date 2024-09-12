import warnings
from typing import Optional, List, Union, Tuple, Dict

import pandas as pd
import sagepy_connector
from .spectrum import ProcessedSpectrum
from .unimod import variable_unimod_mods_to_set, static_unimod_mods_to_set

psc = sagepy_connector.py_scoring
psc_utils = sagepy_connector.py_utility
from .ion_series import IonType
from .mass import Tolerance
from .database import PeptideIx, IndexedDatabase


class ScoreType:
    def __init__(self, name: str):
        name = name.lower()
        names = {"openms", "sage", "happy"}
        assert name in names, f"Invalid score type: {name}, allowed values are: {names}"
        self.__py_ptr = psc.PyScoreType(name)

    @classmethod
    def from_py_ptr(cls, py_ptr: psc.PyScoreType):
        instance = cls.__new__(cls)
        instance.__py_ptr = py_ptr
        return instance

    def get_py_ptr(self) -> psc.PyScoreType:
        return self.__py_ptr

    def __repr__(self):
        return f"ScoreType({self.__py_ptr.to_str()})"


class PeptideSpectrumMatch:
    def __init__(self,
                 spec_idx: str,
                 peptide_idx: int,
                 proteins: List[str],
                 decoy: bool,
                 hyper_score: float,
                 rank: int,
                 mono_mass_observed: Union[None, float],
                 isotope_error: Union[None, int],
                 average_ppm: Union[None, float],
                 delta_next: Union[None, float],
                 delta_best: Union[None, float],
                 matched_peaks: Union[None, int],
                 longest_b: Union[None, int],
                 longest_y: Union[None, int],
                 longest_y_pct: Union[None, float],
                 missed_cleavages: Union[None, int],
                 matched_intensity_pct: Union[None, float],
                 scored_candidates: Union[None, int],
                 poisson: Union[None, float],
                 sequence: Union[None, str],
                 charge: Union[None, int],
                 retention_time_observed: Union[None, float],
                 retention_time_predicted: Union[None, float],
                 inverse_mobility_observed: Union[None, float],
                 inverse_mobility_predicted: Union[None, float],
                 intensity_ms1: Union[None, float],
                 intensity_ms2: Union[None, float],
                 q_value: Union[None, float],
                 collision_energy: Union[None, float],
                 collision_energy_calibrated: Union[None, float],
                 fragments: Union[None, 'Fragments'] = None,
                 re_score: Union[None, float] = None,
                 cosine_similarity: Union[None, float] = None,
                 file_name: Union[None, str] = None,
                 fragment_charges: Union[None, int] = None,
                 fragment_ion_types: Union[None, List[str]] = None,
                 fragment_ordinals: Union[None, List[int]] = None,
                 fragment_intensities: Union[None, List[float]] = None,
                 fragment_mz_calculated: Union[None, List[float]] = None,
                 fragment_mz_experimental: Union[None, List[float]] = None,
                 mz_calibration_ppm: Union[None, float] = None,
                 projected_rt: Union[None, float] = None,
                 beta_score: Union[None, float] = None,
                 posterior_error_prob: Union[None, float] = None,
                 prosit_intensities: Union[None, List[float]] = None,
                 spectral_entropy_similarity: Union[None, float] = None,
                 spectral_correlation_similarity_pearson: Union[None, float] = None,
                 spectral_correlation_similarity_spearman: Union[None, float] = None,
                 spectral_normalized_intensity_difference: Union[None, float] = None,
                 ):
        self.__py_ptr = psc.PyPeptideSpectrumMatch(
            spec_idx, peptide_idx, proteins, decoy, hyper_score, rank, mono_mass_observed,
            isotope_error, average_ppm, delta_next, delta_best, matched_peaks, longest_b, longest_y,
            longest_y_pct, missed_cleavages, matched_intensity_pct, scored_candidates, poisson,
            sequence, charge,
            retention_time_observed, retention_time_predicted, inverse_mobility_observed, inverse_mobility_predicted,
            intensity_ms1, intensity_ms2, q_value, collision_energy, collision_energy_calibrated, fragments.get_py_ptr(),
            re_score, cosine_similarity, file_name,
            fragment_charges, fragment_ion_types, fragment_ordinals, fragment_intensities, fragment_mz_calculated,
            fragment_mz_experimental, mz_calibration_ppm, projected_rt, beta_score, posterior_error_prob, prosit_intensities,
            spectral_entropy_similarity, spectral_correlation_similarity_pearson, spectral_correlation_similarity_spearman,
            spectral_normalized_intensity_difference
        )

    @property
    def spec_idx(self):
        return self.__py_ptr.spec_idx

    @property
    def peptide_idx(self):
        return self.__py_ptr.peptide_idx

    @property
    def proteins(self):
        return self.__py_ptr.proteins

    @property
    def decoy(self):
        return self.__py_ptr.decoy

    @property
    def hyper_score(self):
        return self.__py_ptr.hyper_score

    @hyper_score.setter
    def hyper_score(self, value):
        self.__py_ptr.hyper_score = value

    @property
    def re_score(self):
        return self.__py_ptr.re_score

    @re_score.setter
    def re_score(self, value):
        self.__py_ptr.re_score = value

    @property
    def rank(self):
        return self.__py_ptr.rank

    @property
    def charge(self):
        return self.__py_ptr.charge

    @property
    def mono_mass_observed(self):
        return self.__py_ptr.mono_mass_observed

    @property
    def isotope_error(self):
        return self.__py_ptr.isotope_error

    @property
    def average_ppm(self):
        return self.__py_ptr.average_ppm

    @property
    def delta_next(self):
        return self.__py_ptr.delta_next

    @property
    def delta_best(self):
        return self.__py_ptr.delta_best

    @property
    def matched_peaks(self):
        return self.__py_ptr.matched_peaks

    @property
    def longest_b(self):
        return self.__py_ptr.longest_b

    @property
    def longest_y(self):
        return self.__py_ptr.longest_y

    @property
    def longest_y_pct(self):
        return self.__py_ptr.longest_y_pct

    @property
    def missed_cleavages(self):
        return self.__py_ptr.missed_cleavages

    @property
    def matched_intensity_pct(self):
        return self.__py_ptr.matched_intensity_pct

    @property
    def scored_candidates(self):
        return self.__py_ptr.scored_candidates

    @property
    def poisson(self):
        return self.__py_ptr.poisson

    @property
    def sequence(self):
        return self.__py_ptr.peptide_sequence

    @property
    def mono_mz_calculated(self):
        return self.__py_ptr.mono_mz_calculated

    @property
    def mono_mass_calculated(self):
        return self.__py_ptr.mono_mass_calculated

    @property
    def retention_time_observed(self):
        return self.__py_ptr.retention_time_observed

    @property
    def retention_time_predicted(self):
        return self.__py_ptr.retention_time_predicted

    @retention_time_predicted.setter
    def retention_time_predicted(self, value):
        self.__py_ptr.retention_time_predicted = value

    @property
    def inverse_mobility_observed(self):
        return self.__py_ptr.inverse_mobility_observed

    @property
    def inverse_mobility_predicted(self):
        return self.__py_ptr.inverse_mobility_predicted

    @inverse_mobility_predicted.setter
    def inverse_mobility_predicted(self, value):
        self.__py_ptr.inverse_mobility_predicted = value

    @property
    def intensity_ms1(self):
        return self.__py_ptr.intensity_ms1

    @property
    def intensity_ms2(self):
        return self.__py_ptr.intensity_ms2

    @property
    def q_value(self):
        return self.__py_ptr.q_value

    @property
    def collision_energy(self):
        return self.__py_ptr.collision_energy

    @collision_energy.setter
    def collision_energy(self, value):
        self.__py_ptr.collision_energy = value

    @property
    def collision_energy_calibrated(self):
        return self.__py_ptr.collision_energy_calibrated

    @collision_energy_calibrated.setter
    def collision_energy_calibrated(self, value):
        self.__py_ptr.collision_energy_calibrated = value

    @property
    def cosine_similarity(self):
        return self.__py_ptr.cosine_similarity

    @cosine_similarity.setter
    def cosine_similarity(self, value: float):
        self.__py_ptr.cosine_similarity = value

    @property
    def file_name(self):
        return self.__py_ptr.file_name

    @file_name.setter
    def file_name(self, value):
        self.__py_ptr.file_name = value

    @property
    def fragments_observed(self) -> Union[None, 'Fragments']:
        return Fragments.from_py_fragments(self.__py_ptr.fragments_observed)

    @fragments_observed.setter
    def fragments_observed(self, value: 'Fragments'):
        self.__py_ptr.fragments_observed = value.get_py_ptr()

    @property
    def fragments_predicted(self) -> Union[None, 'Fragments']:
        if self.__py_ptr.fragments_predicted is None:
            return None
        return Fragments.from_py_fragments(self.__py_ptr.fragments_predicted)

    @fragments_predicted.setter
    def fragments_predicted(self, value: 'Fragments'):
        self.__py_ptr.fragments_predicted = value.get_py_ptr()

    @property
    def mz_calibration_ppm(self):
        return self.__py_ptr.mz_calibration_ppm

    @mz_calibration_ppm.setter
    def mz_calibration_ppm(self, value):
        self.__py_ptr.mz_calibration_ppm = value

    @property
    def projected_rt(self):
        return self.__py_ptr.projected_rt

    @projected_rt.setter
    def projected_rt(self, value):
        self.__py_ptr.projected_rt = value

    @property
    def beta_score(self):
        return self.__py_ptr.beta_score

    @beta_score.setter
    def beta_score(self, value):
        self.__py_ptr.beta_score = value

    @property
    def posterior_error_prob(self):
        return self.__py_ptr.posterior_error_prob

    @posterior_error_prob.setter
    def posterior_error_prob(self, value):
        self.__py_ptr.posterior_error_prob = value

    @property
    def prosit_intensities(self):
        return self.__py_ptr.prosit_intensities

    @prosit_intensities.setter
    def prosit_intensities(self, value):
        self.__py_ptr.prosit_intensities = value

    @property
    def spectral_entropy_similarity(self):
        return self.__py_ptr.spectral_entropy_similarity

    @spectral_entropy_similarity.setter
    def spectral_entropy_similarity(self, value):
        self.__py_ptr.spectral_entropy_similarity = value

    @property
    def spectral_correlation_similarity_pearson(self):
        return self.__py_ptr.spectral_correlation_similarity_pearson

    @spectral_correlation_similarity_pearson.setter
    def spectral_correlation_similarity_pearson(self, value):
        self.__py_ptr.spectral_correlation_similarity_pearson = value

    @property
    def spectral_correlation_similarity_spearman(self):
        return self.__py_ptr.spectral_correlation_similarity_spearman

    @spectral_correlation_similarity_spearman.setter
    def spectral_correlation_similarity_spearman(self, value):
        self.__py_ptr.spectral_correlation_similarity_spearman = value

    @property
    def spectral_normalized_intensity_difference(self):
        return self.__py_ptr.spectral_normalized_intensity_difference

    @spectral_normalized_intensity_difference.setter
    def spectral_normalized_intensity_difference(self, value):
        self.__py_ptr.spectral_normalized_intensity_difference = value

    def prosit_fragment_map(self) -> Optional[Dict[Tuple[int, int, int], float]]:
        return self.__py_ptr.prosit_fragment_map()

    def observed_fragment_map(self) -> Optional[Dict[Tuple[int, int, int], float]]:
        return self.__py_ptr.observed_fragment_map()

    @classmethod
    def from_py_ptr(cls, py_ptr: psc.PyPeptideSpectrumMatch):
        instance = cls.__new__(cls)
        instance.__py_ptr = py_ptr
        return instance

    @staticmethod
    def from_json(json_str: str) -> 'PeptideSpectrumMatch':
        psm = psc.psm_from_json(json_str)
        return PeptideSpectrumMatch.from_py_ptr(psm)

    def get_py_ptr(self) -> psc.PyPeptideSpectrumMatch:
        return self.__py_ptr

    def associate_fragment_ions_with_prosit_predicted_intensities(self, flat_intensities: List[float]):
        self.__py_ptr.associate_fragment_ions_with_prosit_predicted_intensities(flat_intensities)

    def match_observed_predicted_intensities(self) -> Tuple[Optional['Fragments'], Optional[List[float]]]:
        maybe_fragment, flat_intensities = self.__py_ptr.match_observed_predicted_intensities()
        if maybe_fragment is not None:
            return Fragments.from_py_fragments(maybe_fragment), flat_intensities
        return None, None

    def to_json(self) -> str:
        return self.__py_ptr.to_json()

    def __repr__(self):
       return (f"PeptideSpectrumMatch(spec_idx: {self.spec_idx}, match_idx: {self.peptide_idx}, "
               f"proteins: {self.proteins}, decoy: {self.decoy}, hyper_score: {self.hyper_score}, "
               f"rank: {self.rank}, mono_mass_observed: {self.mono_mass_observed}, "
               f"isotope_error: {self.isotope_error}, average_ppm: {self.average_ppm}, "
               f"delta_next: {self.delta_next}, delta_best: {self.delta_best}, "
               f"matched_peaks: {self.matched_peaks}, longest_b: {self.longest_b}, "
               f"longest_y: {self.longest_y}, longest_y_pct: {self.longest_y_pct}, "
               f"missed_cleavages: {self.missed_cleavages}, matched_intensity_pct: {self.matched_intensity_pct}, "
               f"scored_candidates: {self.scored_candidates}, poisson: {self.poisson}, "
               f"sequence: {self.sequence}, charge: {self.charge}, "
               f"retention_time_observed: {self.retention_time_observed}, "
               f"retention_time_predicted: {self.retention_time_predicted}, "
               f"inverse_mobility_observed: {self.inverse_mobility_observed}, "
               f"inverse_mobility_predicted: {self.inverse_mobility_predicted}, "
               f"intensity_ms1: {self.intensity_ms1}, intensity_ms2: {self.intensity_ms2}, "
               f"q_value: {self.q_value}, collision_energy: {self.collision_energy}, "
               f"collision_energy_calibrated: {self.collision_energy_calibrated}, "
               f"fragments_observed: {self.fragments_observed if self.fragments_predicted is not None else None}, "
               f"fragments_predicted: {self.fragments_predicted if self.fragments_observed is not None else None}, "
               f"re_score: {self.re_score if self.re_score is not None else None}, "
               f"cosine_similarity: {self.cosine_similarity if self.cosine_similarity is not None else None}, "
               f"file_name: {self.file_name if self.file_name is not None else None}, "
               f"mz_calibration_ppm: {self.mz_calibration_ppm if self.mz_calibration_ppm is not None else None}, "
               f"projected_rt: {self.projected_rt if self.projected_rt is not None else None}, "
               f"beta_score: {self.beta_score if self.beta_score is not None else None}, "
               f"posterior_error_prob: {self.posterior_error_prob if self.posterior_error_prob is not None else None})")


class Fragments:
    def __init__(self,
                charges: List[int],
                ion_types: List[IonType],
                fragment_ordinals: List[int],
                intensities: List[float],
                mz_calculated: List[float],
                mz_experimental: List[float]
                ):

        kinds = [x.get_py_ptr() for x in ion_types]
        self.__fragments_ptr = psc.PyFragments(charges, kinds, fragment_ordinals,
                                               intensities, mz_calculated, mz_experimental)

    @classmethod
    def from_py_fragments(cls, fragments: psc.PyFragments):
        instance = cls.__new__(cls)
        instance.__fragments_ptr = fragments
        return instance

    @property
    def charges(self) -> List[int]:
        return self.__fragments_ptr.charges

    @property
    def ion_types(self) -> List[IonType]:
        return [IonType.from_py_kind(x) for x in self.__fragments_ptr.kinds]

    @property
    def fragment_ordinals(self) -> List[int]:
        return self.__fragments_ptr.fragment_ordinals

    @property
    def intensities(self) -> List[float]:
        return self.__fragments_ptr.intensities

    @property
    def mz_calculated(self) -> List[float]:
        return self.__fragments_ptr.mz_calculated

    @property
    def mz_experimental(self) -> List[float]:
        return self.__fragments_ptr.mz_experimental

    def __repr__(self):
        return (f"Fragments(charges: {self.charges}, "
                f"ion_types: {self.ion_types}, "
                f"fragment_ordinals: {self.fragment_ordinals}, "
                f"intensities: {self.intensities}, "
                f"mz_calculated: {self.mz_calculated}, "
                f"mz_experimental: {self.mz_experimental})")

    def get_py_ptr(self):
        return self.__fragments_ptr


class Scorer:

    def __init__(
            self,
            precursor_tolerance: Tolerance = Tolerance(da=(-5, 5)),
            fragment_tolerance: Tolerance = Tolerance(ppm=(-10, 10)),
            min_matched_peaks: int = 6,
            min_isotope_err: int = -1,
            max_isotope_err: int = 3,
            min_precursor_charge: int = 2,
            max_precursor_charge: int = 4,
            min_fragment_mass: float = 150,
            max_fragment_mass: float = 2000,
            chimera: bool = False,
            report_psms: int = 1,
            wide_window: bool = False,
            annotate_matches: bool = True,
            score_type: ScoreType = ScoreType("openms"),
            max_fragment_charge: Optional[int] = 1,
            variable_mods: Union[Dict[str, List[str]], Dict[str, List[int]]] = None,
            static_mods: Union[Dict[str, str], Dict[str, int]] = None,
    ):
        """Scorer class

        Args:
            precursor_tolerance (Tolerance, optional): The precursor tolerance. Defaults to Tolerance(da=(-5, 5)).
            fragment_tolerance (Tolerance, optional): The fragment tolerance. Defaults to Tolerance(ppm=(-10, 10)).
            min_matched_peaks (int, optional): The minimum number of matched peaks. Defaults to 6.
            min_isotope_err (int, optional): The minimum isotope error. Defaults to -1.
            max_isotope_err (int, optional): The maximum isotope error. Defaults to 3.
            min_precursor_charge (int, optional): The minimum precursor charge. Defaults to 2.
            max_precursor_charge (int, optional): The maximum precursor charge. Defaults to 4.
            min_fragment_mass (float, optional): The minimum fragment mass. Defaults to 150.
            max_fragment_mass (float, optional): The maximum fragment mass. Defaults to 2000.
            chimera (bool, optional): Should chimera be used. Defaults to False.
            report_psms (int, optional): The number of PSMs to report. Defaults to 1.
            wide_window (bool, optional): Should wide window be used. Defaults to False.
            annotate_matches (bool, optional): Should matches be annotated. Defaults to True.
            score_type (ScoreType, optional): The score type. Defaults to ScoreType("openms").
            max_fragment_charge (Optional[int], optional): The maximum fragment charge. Defaults to 1.
        """

        if variable_mods is not None:
            variable_mods = variable_unimod_mods_to_set(variable_mods)
        else:
            warnings.warn("CAUTION! No variable modifications provided, using an empty set.")
            variable_mods = set()

        if static_mods is not None:
            static_mods = static_unimod_mods_to_set(static_mods)
        else:
            warnings.warn("CAUTION! No static modifications provided, using an empty set.")
            static_mods = set()

        expected_modifications = variable_mods.union(static_mods)

        self.__scorer_ptr = psc.PyScorer(
            precursor_tolerance.get_py_ptr(),
            fragment_tolerance.get_py_ptr(),
            min_matched_peaks,
            min_isotope_err,
            max_isotope_err,
            min_precursor_charge,
            max_precursor_charge,
            min_fragment_mass,
            max_fragment_mass,
            chimera,
            report_psms,
            wide_window,
            annotate_matches,
            expected_modifications,
            max_fragment_charge,
            score_type.get_py_ptr()
        )

    @classmethod
    def from_py_scorer(cls, scorer: psc.PyScorer):
        instance = cls.__new__(cls)
        instance.__scorer_ptr = scorer
        return instance

    @property
    def precursor_tolerance(self) -> Tolerance:
        return Tolerance.from_py_tolerance(self.__scorer_ptr.precursor_tolerance)

    @property
    def fragment_tolerance(self) -> Tolerance:
        return Tolerance.from_py_tolerance(self.__scorer_ptr.fragment_tolerance)

    @property
    def min_matched_peaks(self) -> int:
        return self.__scorer_ptr.min_matched_peaks

    @property
    def min_isotope_err(self) -> int:
        return self.__scorer_ptr.min_isotope_err

    @property
    def max_isotope_err(self) -> int:
        return self.__scorer_ptr.max_isotope_err

    @property
    def min_precursor_charge(self) -> int:
        return self.__scorer_ptr.min_precursor_charge

    @property
    def max_precursor_charge(self) -> int:
        return self.__scorer_ptr.max_precursor_charge

    @property
    def min_fragment_mass(self) -> float:
        return self.__scorer_ptr.min_fragment_mass

    @property
    def max_fragment_mass(self) -> float:
        return self.__scorer_ptr.max_fragment_mass

    @property
    def chimera(self) -> bool:
        return self.__scorer_ptr.chimera

    @property
    def report_psms(self) -> int:
        return self.__scorer_ptr.report_psms

    @property
    def wide_window(self) -> bool:
        return self.__scorer_ptr.wide_window

    @property
    def max_fragment_charge(self) -> Optional[int]:
        if self.__scorer_ptr.max_fragment_charge is None:
            return None
        else:
            return self.__scorer_ptr.max_fragment_charge

    def __repr__(self):
        return (f"Scorer({self.precursor_tolerance}, {self.fragment_tolerance}, {self.min_matched_peaks}, "
                f"{self.min_isotope_err}, {self.max_isotope_err}, {self.min_precursor_charge}, "
                f"{self.max_precursor_charge}, {self.min_fragment_mass}, {self.max_fragment_mass}, "
                f"{self.chimera}, {self.report_psms}, {self.wide_window}, {self.max_fragment_charge})")

    def score(self, db: IndexedDatabase, spectrum: ProcessedSpectrum) -> List['Feature']:
        return [Feature.from_py_feature(f) for f in self.__scorer_ptr.score(db.get_py_ptr(), spectrum.get_py_ptr())]

    def score_collection_top_n(self, db: IndexedDatabase,
                               spectrum_collection: List[ProcessedSpectrum], num_threads: int = 4) -> List[
        List['Feature']]:
        scores = self.__scorer_ptr.score_collection(db.get_py_ptr(),
                                                    [spec.get_py_ptr() for spec in spectrum_collection], num_threads)
        return [[Feature.from_py_feature(f) for f in score] for score in scores]

    def score_collection(self, db: IndexedDatabase, spectrum_collection: List[Optional[ProcessedSpectrum]],
                         num_threads: int = 4) -> List['Feature']:
        scores = self.score_collection_top_n(db, spectrum_collection, num_threads)

        result = []

        for score in scores:
            if len(score) > 0:
                result.append(score[0])
            else:
                result.append(None)

        return result

    def score_collection_psm(self, db: IndexedDatabase, spectrum_collection: List[Optional[ProcessedSpectrum]],
                             num_threads: int = 4) -> Dict[str, List[PeptideSpectrumMatch]]:

        py_psms = self.__scorer_ptr.score_collection_to_psm_collection(db.get_py_ptr(),
                                                                       [spec.get_py_ptr() for spec in
                                                                        spectrum_collection],
                                                                       num_threads)
        ret_dict = {}
        for key, values in py_psms.items():
            ret_dict[key] = [PeptideSpectrumMatch.from_py_ptr(psm) for psm in values]

        return ret_dict

    def score_psm(self, db: IndexedDatabase, spectrum: ProcessedSpectrum) -> Dict[str, List[PeptideSpectrumMatch]]:
        py_psms = self.__scorer_ptr.score_collection_to_psm_collection(db.get_py_ptr(), [spectrum.get_py_ptr()], 1)

        ret_dict = {}

        for key, values in py_psms.items():
            ret_dict[key] = [PeptideSpectrumMatch.from_py_ptr(psm) for psm in values]

        return ret_dict

    def _score_chimera_fast(self, db: IndexedDatabase, spectrum: ProcessedSpectrum) -> List['Feature']:
        return [Feature.from_py_feature(f) for f in
                self.__scorer_ptr.score_chimera_fast(db.get_py_ptr(), spectrum.get_py_ptr())]

    def _score_standard(self, db: IndexedDatabase, spectrum: ProcessedSpectrum) -> List['Feature']:
        return [Feature.from_py_feature(f) for f in
                self.__scorer_ptr.score_standard(db.get_py_ptr(), spectrum.get_py_ptr())]


class Feature:
    def __init__(self, peptide_idx: PeptideIx, psm_id: int, peptide_len: int, spec_id: str, file_id: int,
                 rank: int, label: int, expmass: float, calcmass: float, charge: int, delta_mass: float,
                 isotope_error: float, average_ppm: float, hyperscore: float, delta_next: float,
                 delta_best: float, matched_peaks: int, longest_b: int, longest_y: int,
                 longest_y_pct: float, missed_cleavages: int, matched_intensity_pct: float,
                 scored_candidates: int, poisson: float, discriminant_score: float,
                 posterior_error: float, spectrum_q: float, peptide_q: float, protein_q: float,
                 ms2_intensity: float,
                 fragments: Optional[float] = None,
                 rt: Optional[float] = None,
                 aligned_rt: Optional[float] = None,
                 predicted_rt: Optional[float] = None,
                 delta_rt_model: Optional[float] = None,
                 ims: Optional[float] = None,
                 predicted_ims: Optional[float] = None,
                 delta_ims_model: Optional[float] = None, ):
        """Feature class

        Args:
            peptide_idx (PeptideIx): The peptide index
            psm_id (int): The PSM id
            peptide_len (int): The peptide length
            spec_id (str): The spectrum id
            file_id (int): The file id
            rank (int): The rank
            label (int): The label
            expmass (float): The experimental mass
            calcmass (float): The calculated mass
            charge (int): The charge
            rt (float): The retention time
            aligned_rt (float): The aligned retention time
            predicted_rt (float): The predicted retention time
            delta_rt_model (float): The delta retention time model
            delta_mass (float): The delta mass
            isotope_error (float): The isotope error
            average_ppm (float): The average ppm
            hyperscore (float): The hyperscore
            delta_next (float): The delta next
            delta_best (float): The delta best
            matched_peaks (int): The number of matched peaks
            longest_b (int): The longest b ion
            longest_y (int): The longest y ion
            longest_y_pct (float): The longest y ion percentage
            missed_cleavages (int): The number of missed cleavages
            matched_intensity_pct (float): The matched intensity percentage
            scored_candidates (int): The number of scored candidates
            poisson (float): The poisson
            discriminant_score (float): The discriminant score
            posterior_error (float): The posterior error
            spectrum_q (float): The spectrum q
            peptide_q (float): The peptide q
            protein_q (float): The protein q
            ms2_intensity (float): The MS2 intensity
        """

        self.__feature_ptr = psc.PyFeature(peptide_idx, psm_id, peptide_len, spec_id, file_id, rank, label,
                                           expmass, calcmass, charge, delta_mass, isotope_error, average_ppm,
                                           hyperscore, delta_next, delta_best, matched_peaks,
                                           longest_b, longest_y, longest_y_pct, missed_cleavages,
                                           matched_intensity_pct, scored_candidates, poisson,
                                           discriminant_score, posterior_error, spectrum_q,
                                           peptide_q, protein_q, ms2_intensity,
                                           fragments.get_py_ptr() if fragments is not None else None,
                                           rt,
                                           aligned_rt,
                                           predicted_rt,
                                           delta_rt_model,
                                           ims,
                                           predicted_ims,
                                           delta_ims_model, )

    @classmethod
    def from_py_feature(cls, feature: psc.PyFeature):
        instance = cls.__new__(cls)
        instance.__feature_ptr = feature
        return instance

    @property
    def peptide_idx(self) -> PeptideIx:
        return PeptideIx.from_py_peptide_ix(self.__feature_ptr.peptide_idx)

    @property
    def peptide_len(self) -> int:
        return self.__feature_ptr.peptide_len

    @property
    def spec_id(self) -> str:
        return self.__feature_ptr.spec_id

    @property
    def psm_id(self) -> int:
        return self.__feature_ptr.psm_id

    @property
    def file_id(self) -> int:
        return self.__feature_ptr.file_id

    @property
    def rank(self) -> int:
        return self.__feature_ptr.rank

    @property
    def label(self) -> int:
        return self.__feature_ptr.label

    @property
    def expmass(self) -> float:
        return self.__feature_ptr.expmass

    @property
    def calcmass(self) -> float:
        return self.__feature_ptr.calcmass

    @property
    def charge(self) -> int:
        return self.__feature_ptr.charge

    @property
    def rt(self) -> float:
        return self.__feature_ptr.rt

    @property
    def aligned_rt(self) -> float:
        return self.__feature_ptr.aligned_rt

    @property
    def predicted_rt(self) -> float:
        return self.__feature_ptr.predicted_rt

    @property
    def delta_rt_model(self) -> float:
        return self.__feature_ptr.delta_rt_model

    @property
    def delta_mass(self) -> float:
        return self.__feature_ptr.delta_mass

    @property
    def isotope_error(self) -> float:
        return self.__feature_ptr.isotope_error

    @property
    def average_ppm(self) -> float:
        return self.__feature_ptr.average_ppm

    @property
    def hyperscore(self) -> float:
        return self.__feature_ptr.hyperscore

    @property
    def delta_next(self) -> float:
        return self.__feature_ptr.delta_next

    @property
    def delta_best(self) -> float:
        return self.__feature_ptr.delta_best

    @property
    def matched_peaks(self) -> int:
        return self.__feature_ptr.matched_peaks

    @property
    def longest_b(self) -> int:
        return self.__feature_ptr.longest_b

    @property
    def longest_y(self) -> int:
        return self.__feature_ptr.longest_y

    @property
    def longest_y_pct(self) -> float:
        return self.__feature_ptr.longest_y_pct

    @property
    def missed_cleavages(self) -> int:
        return self.__feature_ptr.missed_cleavages

    @property
    def matched_intensity_pct(self) -> float:
        return self.__feature_ptr.matched_intensity_pct

    @property
    def scored_candidates(self) -> int:
        return self.__feature_ptr.scored_candidates

    @property
    def poisson(self) -> float:
        return self.__feature_ptr.poisson

    @property
    def discriminant_score(self) -> float:
        return self.__feature_ptr.discriminant_score

    @property
    def posterior_error(self) -> float:
        return self.__feature_ptr.posterior_error

    @property
    def spectrum_q(self) -> float:
        return self.__feature_ptr.spectrum_q

    @property
    def peptide_q(self) -> float:
        return self.__feature_ptr.peptide_q

    @property
    def protein_q(self) -> float:
        return self.__feature_ptr.protein_q

    @property
    def ms2_intensity(self) -> float:
        return self.__feature_ptr.ms2_intensity

    @property
    def fragments(self) -> Optional[Fragments]:
        if self.__feature_ptr.fragments is None:
            return None
        else:
            return Fragments.from_py_fragments(self.__feature_ptr.fragments)

    @property
    def ims(self) -> Optional[float]:
        return self.__feature_ptr.ims

    @property
    def predicted_ims(self) -> Optional[float]:
        return self.__feature_ptr.predicted_ims

    @property
    def delta_ims_model(self) -> Optional[float]:
        return self.__feature_ptr.delta_ims_model

    def __repr__(self):
        return (f"Feature("
                f"idx: {self.peptide_idx}, "
                f"psm_id: {self.psm_id}, "
                f"peptide_len: {self.peptide_len}, "
                f"spec_id: {self.spec_id}, "
                f"file_id: {self.file_id}, "
                f"rank: {self.rank}, "
                f"label: {self.label}, "
                f"exp. mass: {self.expmass}, "
                f"cal. mass: {self.calcmass}, "
                f"charge: {self.charge}, "
                f"retention time: {self.rt}, "
                f"aligned rt: {self.aligned_rt}, "
                f"predicted rt: {self.predicted_rt}, "
                f"delta rt model: {self.delta_rt_model}, "
                f"delta mass: {self.delta_mass}, "
                f"isotope error: {self.isotope_error}, "
                f"average ppm: {self.average_ppm}, "
                f"hyperscore: {self.hyperscore}, "
                f"delta_next: {self.delta_next}, "
                f"delta_best: {self.delta_best}, "
                f"matched peaks: {self.matched_peaks}, "
                f"longest b: {self.longest_b},"
                f"longest y: {self.longest_y}, "
                f"longest y pct: {self.longest_y_pct}, "
                f"missed cleavages: {self.missed_cleavages}, "
                f"matched intensity pct: {self.matched_intensity_pct}, "
                f"scored candidates: {self.scored_candidates}, "
                f"poisson: {self.poisson}, "
                f"discriminant score: {self.discriminant_score}, "
                f"posterior error: {self.posterior_error}, "
                f"spectrum q: {self.spectrum_q}, "
                f"peptide q: {self.peptide_q}, "
                f"protein q: {self.protein_q}, "
                f"ms2 intensity: {self.ms2_intensity}), "
                f"fragments: {self.fragments})")

    def get_py_ptr(self):
        return self.__feature_ptr


def associate_fragment_ions_with_prosit_predicted_intensities(
        psms: List[PeptideSpectrumMatch],
        flat_intensities: List[List[float]], num_threads: int = 16) -> List['PeptideSpectrumMatch']:
    """Associate fragment ions with prosit predicted intensities in parallel

    Args:
        psms (List[PeptideSpectrumMatch]): The peptide spectrum matches
        flat_intensities (List[List[float]]): The flat intensities
        num_threads (int, optional): The number of threads. Defaults to 16.

    Returns:
        List[PeptideSpectrumMatch]: The peptide spectrum matches
    """
    result = psc.associate_fragment_ions_with_prosit_predicted_intensities_par(
        [psm.get_py_ptr() for psm in psms], flat_intensities, num_threads
    )
    return [PeptideSpectrumMatch.from_py_ptr(f) for f in result]


def associate_fragment_ions_with_prosit_predicted_intensities_pandas(
        psms: List[PeptideSpectrumMatch],
        flat_intensities: List[List[float]],
        num_threads: int = 16,
) -> pd.DataFrame:
    """Associate fragment ions with prosit predicted intensities in parallel

    Args:
        psms (List[PeptideSpectrumMatch]): The peptide spectrum matches
        flat_intensities (List[List[float]]): The flat intensities
        num_threads (int, optional): The number of threads. Defaults to 16.

    Returns:
        pd.DataFrame: The pandas dataframe
    """
    result = psc.associate_fragment_ions_with_prosit_predicted_intensities_par(
        [psm.get_py_ptr() for psm in psms], flat_intensities, num_threads
    )
    row_list = []
    for match in result:
        row_list.append({
            "spec_idx": match.spec_idx,
            "match_idx": match.peptide_idx,
            "proteins": match.proteins,
            "decoy": match.decoy,
            "score": match.hyper_score,
            "rank": match.rank,
            "mono_mz_calculated": match.mono_mz_calculated,
            "mono_mass_observed": match.mono_mass_observed,
            "mono_mass_calculated": match.mono_mass_calculated,
            "sequence": match.peptide_sequence,
            "charge": match.charge,
            "retention_time_observed": match.retention_time_observed,
            "retention_time_predicted": match.retention_time_predicted,
            "inverse_mobility_observed": match.inverse_mobility_observed,
            "inverse_mobility_predicted": match.inverse_mobility_predicted,
            "intensity_ms1": match.intensity_ms1,
            "intensity_ms2": match.intensity_ms2,
            "q_value": match.q_value,
            "collision_energy": match.collision_energy,
            "cosine_similarity": match.cosine_similarity,
        })
    return pd.DataFrame(row_list)


def json_bin_to_psms(json_bin: bytes) -> List[PeptideSpectrumMatch]:
    """ Convert a binary JSON string to a list of PeptideSpectrumMatch objects.

    Args:
        json_bin: a binary JSON string

    Returns:
        a list of PeptideSpectrumMatch objects
    """
    return [PeptideSpectrumMatch.from_py_ptr(json_str) for json_str in psc_utils.json_bin_to_psms(json_bin)]


def psms_to_json(psms, num_threads: int = 4) -> List[str]:
    """ Convert a list of PeptideSpectrumMatch objects to a JSON string.

    Args:
        psms: a list of PeptideSpectrumMatch objects
        num_threads: the number of threads to use

    Returns:
        a JSON string
    """
    return psc_utils.psms_to_json([psm.get_py_ptr() for psm in psms], num_threads)


def psms_to_json_bin(psms) -> bytes:
    """ Convert a list of PeptideSpectrumMatch objects to a binary JSON string.

    Args:
        psms: a list of PeptideSpectrumMatch objects

    Returns:
        a binary JSON string
    """
    return psc_utils.psms_to_json_bin([psm.get_py_ptr() for psm in psms])


def merge_psm_dicts(left_psms: Dict[str, List[PeptideSpectrumMatch]],
                    right_psms: Dict[str, List[PeptideSpectrumMatch]],
                    max_hits: int = 5) -> Dict[str, List[PeptideSpectrumMatch]]:
    """ Merge two dictionaries of peptide spectrum matches.

    Args:
        left_psms: the left dictionary of peptide spectrum matches
        right_psms: the right dictionary of peptide spectrum matches
        max_hits: the maximum number of hits

    Returns:
        a dictionary of peptide spectrum matches
    """
    left_map = {key: [psm.get_py_ptr() for psm in value] for key, value in left_psms.items()}
    right_map = {key: [psm.get_py_ptr() for psm in value] for key, value in right_psms.items()}
    result = psc.merge_psm_maps(left_map, right_map, max_hits)
    return {key: [PeptideSpectrumMatch.from_py_ptr(psm) for psm in value] for key, value in result.items()}


def prosit_intensities_to_fragments(
        flat_intensities: List[float],
) -> Fragments:
    """ Convert a list of intensities to a Fragments object.

    Args:
        flat_intensities: a list of intensities

    Returns:
        a Fragments object
    """
    return Fragments.from_py_fragments(psc.prosit_intensities_to_py_fragments(flat_intensities))


def prosit_intensities_to_fragments_par(
        flat_intensities: List[float],
        num_threads: int = 16,
) -> List[Fragments]:
    """ Convert a list of intensities to a Fragments object in parallel.

    Args:
        flat_intensities: a list of intensities
        num_threads: the number of threads

    Returns:
        a Fragments List object
    """

    return [Fragments.from_py_fragments(f) for f in psc.prosit_intensities_to_py_fragments_par(flat_intensities, num_threads)]
