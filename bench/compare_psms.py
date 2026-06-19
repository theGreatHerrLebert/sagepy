"""Compare sagepy vs sage CLI PSM tables on the same .d input.

Spectrum-level join is unreliable (sagepy's `spec_idx` indexes the processed
list, sage CLI's `scannr` indexes the raw timsrust reader), so this script
reports peptide-level identification overlap and scoring distributions instead.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent


def load_sagepy(path: Path) -> pd.DataFrame:
    df = pd.read_parquet(path)
    return df.rename(
        columns={
            "decoy": "is_decoy",
            "sequence_modified": "peptide",
            "sequence": "stripped_peptide",
        }
    )


def load_sage(path: Path) -> pd.DataFrame:
    return pd.read_parquet(path)


def filter_top_targets_at_q(df: pd.DataFrame, q_max: float = 0.01) -> pd.DataFrame:
    return df[
        (df["rank"] == 1)
        & (df["is_decoy"] == False)  # noqa: E712 — pandas bool mask
        & (df["spectrum_q"] <= q_max)
    ].copy()


def load_sagepy_timings(path: Path) -> dict | None:
    if not path.exists():
        return None
    with open(path) as f:
        return json.load(f)


def load_sage_time_log(path: Path) -> dict | None:
    if not path.exists():
        return None
    text = path.read_text()
    out: dict[str, float | int | str] = {}

    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("Elapsed (wall clock) time"):
            out["wall_clock"] = stripped.rsplit(": ", 1)[1]

    if match := re.search(r"finished in\s+([0-9.]+)s", text):
        out["reported_total_seconds"] = float(match.group(1))
    if match := re.search(r"Maximum resident set size \(kbytes\):\s*([0-9]+)", text):
        out["max_rss_kb"] = int(match.group(1))
    if match := re.search(r"- file IO:\s+([0-9]+) ms", text):
        out["read_seconds"] = int(match.group(1)) / 1000.0
    if match := re.search(r"- search:\s+([0-9]+) ms", text):
        out["search_seconds"] = int(match.group(1)) / 1000.0
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sagepy", type=Path, default=HERE / "results" / "sagepy_psms.parquet")
    ap.add_argument("--sage", type=Path, default=HERE / "results" / "sage" / "results.sage.parquet")
    ap.add_argument("--sagepy-timings", type=Path, default=HERE / "results" / "sagepy_timings.json")
    ap.add_argument("--sage-log", type=Path, default=HERE / "results" / "sage" / "sage_run.log")
    ap.add_argument("--q", type=float, default=0.01)
    args = ap.parse_args()

    sp = load_sagepy(args.sagepy)
    sg = load_sage(args.sage)

    print(f"raw rows  sagepy={len(sp):>7}  sage={len(sg):>7}")
    print(f"top-1     sagepy={(sp['rank']==1).sum():>7}  sage={(sg['rank']==1).sum():>7}")
    print(f"targets   sagepy={(~sp['is_decoy'].astype(bool)).sum():>7}  sage={(~sg['is_decoy'].astype(bool)).sum():>7}")

    sp_top = filter_top_targets_at_q(sp, args.q)
    sg_top = filter_top_targets_at_q(sg, args.q)
    print(f"top1 target q<{args.q}  sagepy={len(sp_top):>7}  sage={len(sg_top):>7}")

    sagepy_timings = load_sagepy_timings(args.sagepy_timings)
    sage_timings = load_sage_time_log(args.sage_log)
    if sagepy_timings or sage_timings:
        print()
        print("speed and memory")
        if sagepy_timings:
            phase = sagepy_timings["phase_seconds"]
            print(
                "  sagepy "
                f"total={phase['total']:.2f}s "
                f"read={phase['read']:.2f}s "
                f"digest={phase['digest']:.2f}s "
                f"score={phase['score']:.2f}s "
                f"rss={sagepy_timings['max_rss_kb'] / 1024 / 1024:.2f} GiB"
            )
        if sage_timings:
            wall = sage_timings.get("wall_clock", "n/a")
            total = sage_timings.get("reported_total_seconds", "n/a")
            rss = sage_timings.get("max_rss_kb")
            rss_text = f"{rss / 1024 / 1024:.2f} GiB" if isinstance(rss, int) else "n/a"
            print(f"  sage   total={total}s wall={wall} rss={rss_text}")

    # Peptide-level overlap (sequence_modified)
    sp_peps = set(sp_top["peptide"])
    sg_peps = set(sg_top["peptide"])
    inter = sp_peps & sg_peps
    print()
    print("unique peptides at 1% FDR (sequence_modified)")
    print(f"  sagepy: {len(sp_peps):>6}")
    print(f"  sage  : {len(sg_peps):>6}")
    print(f"  shared: {len(inter):>6}")
    print(f"  sagepy-only: {len(sp_peps - sg_peps):>6}")
    print(f"  sage-only  : {len(sg_peps - sp_peps):>6}")

    # Stripped-peptide overlap (treat mod variants as the same peptide)
    sp_strip = set(sp_top["stripped_peptide"])
    sg_strip = set(sg_top["stripped_peptide"])
    inter_strip = sp_strip & sg_strip
    print()
    print("unique peptides at 1% FDR (stripped, mod-agnostic)")
    print(f"  sagepy: {len(sp_strip):>6}")
    print(f"  sage  : {len(sg_strip):>6}")
    print(f"  shared: {len(inter_strip):>6}")
    print(f"  sagepy-only: {len(sp_strip - sg_strip):>6}")
    print(f"  sage-only  : {len(sg_strip - sp_strip):>6}")

    # Hyperscore distribution comparison on the intersection
    print()
    print("hyperscore distribution (top-1 target, q<1%)")
    print(f"  sagepy  mean={sp_top['hyperscore'].mean():.2f}  median={sp_top['hyperscore'].median():.2f}  min={sp_top['hyperscore'].min():.2f}  max={sp_top['hyperscore'].max():.2f}")
    print(f"  sage    mean={sg_top['hyperscore'].mean():.2f}  median={sg_top['hyperscore'].median():.2f}  min={sg_top['hyperscore'].min():.2f}  max={sg_top['hyperscore'].max():.2f}")

    # Sample of sage-only peptides (sequence_modified) to inspect
    sample = sorted(sg_peps - sp_peps)[:15]
    if sample:
        print()
        print("sample sage-only peptides (top 15):")
        for p in sample:
            row = sg_top[sg_top["peptide"] == p].iloc[0]
            print(f"  {p:<40}  scannr={row['scannr']:>6}  charge={int(row['charge'])}  hyperscore={row['hyperscore']:.2f}  matched_peaks={int(row['matched_peaks'])}")


if __name__ == "__main__":
    main()
