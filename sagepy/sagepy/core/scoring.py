from typing import Union, Optional, List

import numpy as np
import sagepy_connector

from .spectrum import ProcessedSpectrum

psc = sagepy_connector.py_scoring

from .mass import Tolerance
from .database import PeptideIx, IndexedDatabase


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
            max_fragment_charge: Optional[int] = 1):
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
            max_fragment_charge (Optional[int], optional): The maximum fragment charge. Defaults to 1.
        """
        self.__scorer_ptr = psc.PyScorer(precursor_tolerance.get_py_ptr(),
                                         fragment_tolerance.get_py_ptr(),
                                         min_matched_peaks,
                                         min_isotope_err, max_isotope_err, min_precursor_charge,
                                         max_precursor_charge, min_fragment_mass, max_fragment_mass,
                                         chimera, report_psms, wide_window, max_fragment_charge)

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
                               spectrum_collection: List[ProcessedSpectrum], num_threads: int = 4) -> List[List['Feature']]:
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

    def _score_chimera_fast(self, db: IndexedDatabase, spectrum: ProcessedSpectrum) -> List['Feature']:
        return [Feature.from_py_feature(f) for f in
                self.__scorer_ptr.score_chimera_fast(db.get_py_ptr(), spectrum.get_py_ptr())]

    def _score_standard(self, db: IndexedDatabase, spectrum: ProcessedSpectrum) -> List['Feature']:
        return [Feature.from_py_feature(f) for f in
                self.__scorer_ptr.score_standard(db.get_py_ptr(), spectrum.get_py_ptr())]


class Feature:
    def __init__(self, peptide_idx: PeptideIx, peptide_len: int, spec_id: str, file_id: int,
                 rank: int, label: int, expmass: float, calcmass: float, charge: int, rt: float,
                 aligned_rt: float, predicted_rt: float, delta_rt_model: float, delta_mass: float,
                 isotope_error: float, average_ppm: float, hyperscore: float, delta_next: float,
                 delta_best: float, matched_peaks: int, longest_b: int, longest_y: int,
                 longest_y_pct: float, missed_cleavages: int, matched_intensity_pct: float,
                 scored_candidates: int, poisson: float, discriminant_score: float,
                 posterior_error: float, spectrum_q: float, peptide_q: float, protein_q: float,
                 ms2_intensity: float, ms1_intensity: float):
        """Feature class

        Args:
            peptide_idx (PeptideIx): The peptide index
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
            ms1_intensity (float): The MS1 intensity
        """

        self.__feature_ptr = psc.PyFeature(peptide_idx, peptide_len, spec_id, file_id, rank, label,
                                           expmass, calcmass, charge, rt, aligned_rt, predicted_rt,
                                           delta_rt_model, delta_mass, isotope_error, average_ppm,
                                           hyperscore, delta_next, delta_best, matched_peaks,
                                           longest_b, longest_y, longest_y_pct, missed_cleavages,
                                           matched_intensity_pct, scored_candidates, poisson,
                                           discriminant_score, posterior_error, spectrum_q,
                                           peptide_q, protein_q, ms2_intensity, ms1_intensity)

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
    def ms1_intensity(self) -> float:
        return self.__feature_ptr.ms1_intensity

    def __repr__(self):
        return (f"Feature("
                f"idx: {self.peptide_idx}, "
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
                f"ms2 intensity: {self.ms2_intensity}, "
                f"ms1 intensity: {self.ms1_intensity})")

    def get_py_ptr(self):
        return self.__feature_ptr
