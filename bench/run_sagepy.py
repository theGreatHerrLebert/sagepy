"""Run a sagepy DDA search on a Bruker .d file using the native TDF reader.

Mirrors the layout of a sage CLI config (bench/sage_config.json) so results
can be compared against the sage CLI directly.

Usage:
    python run_sagepy.py [--config bench/sage_config.json] [--input .../something.d]
"""
from __future__ import annotations

import argparse
import gc
import json
import multiprocessing as mp
import resource
import time
from pathlib import Path
from typing import Iterable

from sagepy.core import (
    EnzymeBuilder,
    SageSearchConfiguration,
    Scorer,
    SpectrumProcessor,
    Tolerance,
)
from sagepy.core.ion_series import IonType
from sagepy.core.fdr import (
    assign_initial_spectrum_q_psms,
    lda_rescore_psms,
    sage_fdr_psm,
)
from sagepy.core.ml.mobility_model import predict_sage_im
from sagepy.core.ml.retention_alignment import global_alignment_psm
from sagepy.core.ml.retention_model import predict_sage_rt
from sagepy.core.scoring import ScoreType
from sagepy.core.spectrum import read_spectra
from sagepy.core.spectrum import read_processed_spectra
from sagepy.core.unimod import unimod_to_mass
from sagepy.utility import psm_collection_to_pandas


HERE = Path(__file__).resolve().parent


def load_config(path: Path) -> dict:
    with open(path) as f:
        return json.load(f)


def build_database(cfg: dict, fasta_path: Path) -> "IndexedDatabase":
    db_cfg = cfg["database"]
    enz = db_cfg["enzyme"]
    enzyme_builder = EnzymeBuilder(
        missed_cleavages=enz["missed_cleavages"],
        min_len=enz["min_len"],
        max_len=enz["max_len"],
        cleave_at=enz["cleave_at"],
        restrict=enz.get("restrict") or "",
        c_terminal=enz.get("c_terminal") if enz.get("c_terminal") is not None else True,
        semi_enzymatic=enz.get("semi_enzymatic"),
    )

    static_mods_raw = db_cfg.get("static_mods", {}) or {}
    variable_mods_raw = db_cfg.get("variable_mods", {}) or {}
    static_mods = {aa: _mass_to_unimod(mass) for aa, mass in static_mods_raw.items()}
    variable_mods = {
        aa: [_mass_to_unimod(m) for m in masses]
        for aa, masses in variable_mods_raw.items()
    }

    with open(fasta_path) as f:
        fasta = f.read()

    ion_kinds_raw = db_cfg.get("ion_kinds")
    ion_kinds = [IonType(k) for k in ion_kinds_raw] if ion_kinds_raw else None

    kwargs = dict(
        fasta=fasta,
        static_mods=static_mods,
        variable_mods=variable_mods,
        enzyme_builder=enzyme_builder,
        generate_decoys=db_cfg.get("generate_decoys", True),
        bucket_size=db_cfg.get("bucket_size", 16384),
        peptide_min_mass=db_cfg.get("peptide_min_mass", 500.0),
        peptide_max_mass=db_cfg.get("peptide_max_mass", 5000.0),
        min_ion_index=db_cfg.get("min_ion_index", 2),
        max_variable_mods=db_cfg.get("max_variable_mods", 2),
        decoy_tag=db_cfg.get("decoy_tag", "rev_"),
    )
    if ion_kinds is not None:
        kwargs["ion_kinds"] = ion_kinds
    for opt in ("prefilter", "prefilter_chunk_size", "prefilter_low_memory"):
        if opt in db_cfg:
            kwargs[opt] = db_cfg[opt]

    sage_config = SageSearchConfiguration(**kwargs)
    return sage_config.generate_indexed_database()


def _build_mass_to_unimod() -> dict[int, str]:
    """Invert sagepy's UNIMOD-to-mass table to mass(rounded)-to-UNIMOD-string."""
    out: dict[int, str] = {}
    for unimod_id, mass in unimod_to_mass().items():
        # 4 decimal places is about 0.1 ppm at 1000 Da, finer than sage configs write.
        out[round(float(mass), 4)] = unimod_id
    return out


_MASS_TO_UNIMOD = _build_mass_to_unimod()


def _mass_to_unimod(value):
    """Translate a sage-style mass (or list) into the UNIMOD string sagepy expects."""
    if isinstance(value, list):
        return [_mass_to_unimod(v) for v in value]
    if isinstance(value, str):
        return value
    key = round(float(value), 4)
    if key in _MASS_TO_UNIMOD:
        return _MASS_TO_UNIMOD[key]
    raise ValueError(
        f"Mass {value!r} has no UNIMOD entry in sagepy's table; check the config or add the mod."
    )


def _ppm_tolerance(cfg_pair) -> Tolerance:
    lo, hi = cfg_pair["ppm"]
    return Tolerance(ppm=(float(lo), float(hi)))


_SCORE_TYPE_MAP = {
    "SageHyperScore": "hyperscore",
    "OpenMSHyperScore": "openmshyperscore",
}


def _score_type(cfg: dict) -> ScoreType:
    raw = cfg.get("score_type", "SageHyperScore")
    return ScoreType(_SCORE_TYPE_MAP.get(raw, raw))


def build_scorer(cfg: dict) -> Scorer:
    static_mods = cfg["database"].get("static_mods", {}) or {}
    variable_mods = cfg["database"].get("variable_mods", {}) or {}
    static_mods = {aa: _mass_to_unimod(mass) for aa, mass in static_mods.items()}
    variable_mods = {
        aa: [_mass_to_unimod(m) for m in masses]
        for aa, masses in variable_mods.items()
    }

    return Scorer(
        precursor_tolerance=_ppm_tolerance(cfg["precursor_tol"]),
        fragment_tolerance=_ppm_tolerance(cfg["fragment_tol"]),
        min_matched_peaks=cfg.get("min_matched_peaks", 4),
        min_isotope_err=cfg["isotope_errors"][0],
        max_isotope_err=cfg["isotope_errors"][1],
        min_precursor_charge=cfg["precursor_charge"][0],
        max_precursor_charge=cfg["precursor_charge"][1],
        chimera=cfg.get("chimera", False),
        report_psms=cfg.get("report_psms", 1),
        wide_window=cfg.get("wide_window", False),
        override_precursor_charge=cfg.get("override_precursor_charge", False),
        score_type=_score_type(cfg),
        max_fragment_charge=cfg.get("max_fragment_charge", 1),
        static_mods=static_mods,
        variable_mods=variable_mods,
    )


def process_spectra(raw_spectra: Iterable, cfg: dict) -> list:
    processor = SpectrumProcessor(
        take_top_n=cfg.get("max_peaks", 150),
        deisotope=cfg.get("deisotope", True),
    )
    min_peaks = int(cfg.get("min_peaks", 15))
    processed = []
    n_dropped_no_precursor = 0
    n_dropped_low_peaks = 0
    for spec in raw_spectra:
        if spec.ms_level != 2 or not spec.precursors:
            n_dropped_no_precursor += 1
            continue
        proc = processor.process(spec)
        if len(proc.peaks) < min_peaks:
            n_dropped_low_peaks += 1
            continue
        processed.append(proc)
    print(f"[process] dropped {n_dropped_no_precursor} (no precursor / non-MS2), "
          f"{n_dropped_low_peaks} (<{min_peaks} peaks)")
    return processed


def process_spectra_parallel(raw_spectra: list, cfg: dict, num_threads: int) -> list:
    min_peaks = int(cfg.get("min_peaks", 15))
    candidates = []
    n_dropped_no_precursor = 0
    for spec in raw_spectra:
        if spec.ms_level != 2 or not spec.precursors:
            n_dropped_no_precursor += 1
        else:
            candidates.append(spec)

    processor = SpectrumProcessor(
        take_top_n=cfg.get("max_peaks", 150),
        deisotope=cfg.get("deisotope", True),
    )
    processed = processor.process_collection(candidates, num_threads=num_threads)
    keep = [spec for spec in processed if len(spec.peaks) >= min_peaks]
    n_dropped_low_peaks = len(processed) - len(keep)
    print(f"[process] dropped {n_dropped_no_precursor} (no precursor / non-MS2), "
          f"{n_dropped_low_peaks} (<{min_peaks} peaks)")
    return keep


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=HERE / "sage_config.json")
    parser.add_argument(
        "--input",
        type=Path,
        default=None,
        help="Override the .d/.mzML/.mgf/.pmsms path (default: first entry in config.mzml_paths)",
    )
    parser.add_argument("--fasta", type=Path, default=HERE / "data" / "human.fasta")
    parser.add_argument(
        "--output",
        type=Path,
        default=HERE / "results" / "sagepy_psms.parquet",
    )
    parser.add_argument("--num-threads", type=int, default=mp.cpu_count())
    parser.add_argument(
        "--requires-ms1",
        action="store_true",
        help="Also load MS1 frames when reading .d (slower, off by default).",
    )
    parser.add_argument(
        "--rescore",
        action="store_true",
        help="Fit RT + IM models, then run FDR on the sage discriminant score "
        "(matches sage CLI's predict_rt=true behavior).",
    )
    parser.add_argument(
        "--fused-processing",
        action="store_true",
        help="Read and process MS2 spectra in native Rust before returning to Python.",
    )
    parser.add_argument(
        "--timings-output",
        type=Path,
        default=HERE / "results" / "sagepy_timings.json",
        help="Write phase timings and summary counters as JSON.",
    )
    args = parser.parse_args()

    cfg = load_config(args.config)
    input_path = args.input or Path(cfg["mzml_paths"][0])
    args.output.parent.mkdir(parents=True, exist_ok=True)

    timings: dict[str, float] = {}

    print(f"[read]    input: {input_path}")
    if args.fused_processing:
        t0 = time.perf_counter()
        processed_all = read_processed_spectra(
            str(input_path),
            take_top_n=cfg.get("max_peaks", 150),
            deisotope=cfg.get("deisotope", True),
            num_threads=args.num_threads,
            requires_ms1=args.requires_ms1,
        )
        min_peaks = int(cfg.get("min_peaks", 15))
        processed = [spec for spec in processed_all if len(spec.peaks) >= min_peaks]
        n_dropped_low_peaks = len(processed_all) - len(processed)
        del processed_all
        gc.collect()
        timings["process"] = time.perf_counter() - t0
        n_raw_spectra = len(processed) + n_dropped_low_peaks
        n_ms2 = n_raw_spectra
        n_ms1 = 0
        timings["read"] = 0.0
        print(f"[process] dropped 0 (no precursor / non-MS2), "
              f"{n_dropped_low_peaks} (<{min_peaks} peaks)")
        print(f"[process] {len(processed)} processed MS2 in {timings['process']:.2f}s")
    else:
        t0 = time.perf_counter()
        raw = read_spectra(str(input_path), requires_ms1=args.requires_ms1)
        timings["read"] = time.perf_counter() - t0
        n_ms2 = sum(1 for s in raw if s.ms_level == 2)
        n_ms1 = sum(1 for s in raw if s.ms_level == 1)
        print(f"[read]    {len(raw)} spectra ({n_ms2} MS2, {n_ms1} MS1) in {timings['read']:.2f}s")

        t0 = time.perf_counter()
        processed = process_spectra_parallel(raw, cfg, args.num_threads)
        n_raw_spectra = len(raw)
        del raw
        gc.collect()
        timings["process"] = time.perf_counter() - t0
        print(f"[process] {len(processed)} processed MS2 in {timings['process']:.2f}s")

    if not processed:
        raise SystemExit("No MS2 spectra to score; check the input file or filters.")

    t0 = time.perf_counter()
    db = build_database(cfg, args.fasta)
    timings["digest"] = time.perf_counter() - t0
    print(f"[digest]  indexed db built in {timings['digest']:.2f}s")

    scorer = build_scorer(cfg)
    t0 = time.perf_counter()
    psms_dict = scorer.score_collection_psm(
        db=db,
        spectrum_collection=processed,
        num_threads=args.num_threads,
    )
    timings["score"] = time.perf_counter() - t0
    psm_list = [psm for psms in psms_dict.values() for psm in psms]
    print(f"[score]   {len(psm_list)} PSMs across {len(psms_dict)} spectra in {timings['score']:.2f}s")

    if args.rescore:
        t0 = time.perf_counter()
        assign_initial_spectrum_q_psms(psm_list)
        timings["prefdr"] = time.perf_counter() - t0
        print(f"[prefdr]  initial spectrum-q pass in {timings['prefdr']:.2f}s")

        t0 = time.perf_counter()
        alignments = global_alignment_psm(psm_list)
        timings["align"] = time.perf_counter() - t0
        print(f"[align]   global RT alignment over {len(alignments)} file(s) in {timings['align']:.2f}s")

        t0 = time.perf_counter()
        predict_sage_rt(psm_list, db)
        predict_sage_im(psm_list, db)
        timings["rescore"] = time.perf_counter() - t0
        print(f"[rescore] RT + IM models fit in {timings['rescore']:.2f}s")

        t0 = time.perf_counter()
        fitted = lda_rescore_psms(psm_list, _ppm_tolerance(cfg["precursor_tol"]))
        timings["lda"] = time.perf_counter() - t0
        print(f"[lda]     fitted={fitted} in {timings['lda']:.2f}s")

    t0 = time.perf_counter()
    sage_fdr_psm(psm_list, db, use_hyper_score=not args.rescore)
    timings["fdr"] = time.perf_counter() - t0
    print(f"[fdr]     ({'discriminant' if args.rescore else 'hyperscore'}) "
          f"{timings['fdr']:.2f}s")

    t0 = time.perf_counter()
    df = psm_collection_to_pandas(psm_list)
    df.to_parquet(args.output)
    timings["write"] = time.perf_counter() - t0
    print(f"[write]   {len(df)} rows -> {args.output} in {timings['write']:.2f}s")

    timings["total"] = sum(timings.values())
    print()
    print("phase    seconds")
    for phase, secs in timings.items():
        print(f"  {phase:<7} {secs:7.2f}")

    if "spectrum_q" in df.columns:
        below_1pct = (df["spectrum_q"] < 0.01).sum()
        print(f"\nPSMs at 1% spectrum-q FDR: {below_1pct}")

    usage = resource.getrusage(resource.RUSAGE_SELF)
    summary = {
        "input": str(input_path),
        "fasta": str(args.fasta),
        "output": str(args.output),
        "num_threads": args.num_threads,
        "requires_ms1": args.requires_ms1,
        "rescore": args.rescore,
        "n_raw_spectra": n_raw_spectra,
        "n_ms2_raw": n_ms2,
        "n_ms1_raw": n_ms1,
        "n_processed_ms2": len(processed),
        "n_psms": len(df),
        "n_spectra_with_psms": len(psms_dict),
        "psms_spectrum_q_lt_0_01": int((df["spectrum_q"] < 0.01).sum())
        if "spectrum_q" in df.columns
        else None,
        "phase_seconds": timings,
        "max_rss_kb": usage.ru_maxrss,
    }
    args.timings_output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.timings_output, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[timing] wrote {args.timings_output}")


if __name__ == "__main__":
    main()
