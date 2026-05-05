"""Convert Matteo's mmappet+MGF layout into sagepy's .pmsms wrapper.

The full F9477 archive contains:

    optimal2tier/
      pmsms.mmappet/
        0.bin, 1.bin, 2.bin
        dataindex.mmappet/
      plain/mgf.mgf

sagepy's native pmsms reader expects:

    output.pmsms/
      pmsms.mmappet/
      tof2mz.mmappet/
      precursors.parquet

This script derives the missing metadata by streaming the MGF in the same
spectrum order as pmsms.mmappet/dataindex.mmappet. The large pmsms.mmappet
directory is symlinked by default, so the conversion writes only the smaller
tof2mz and precursor metadata.
"""
from __future__ import annotations

import argparse
import os
import re
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd


TITLE_PRECURSOR_RE = re.compile(r"\bprecursor_idx=(\d+)\b")
TITLE_IIM_RE = re.compile(r"\biim=([0-9.eE+-]+)\b")
TITLE_CHARGE_RE = re.compile(r"\bcharge=([0-9]+)\b")


def read_mmappet_schema(path: Path) -> list[tuple[str, str]]:
    schema = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            dtype, name = line.split(maxsplit=1)
            schema.append((dtype, name))
    return schema


def load_dataindex(dataindex_dir: Path) -> dict[str, np.memmap]:
    expected = [
        ("uint64", "precursor_idx", np.uint64),
        ("uint64", "size", np.uint64),
        ("uint64", "idx", np.uint64),
        ("uint32", "max_group_len", np.uint32),
        ("float32", "avg_group_len", np.float32),
    ]
    schema = read_mmappet_schema(dataindex_dir / "schema.txt")
    if schema != [(dtype, name) for dtype, name, _ in expected]:
        raise ValueError(f"Unsupported dataindex schema in {dataindex_dir}: {schema!r}")

    return {
        name: np.memmap(dataindex_dir / f"{i}.bin", mode="r", dtype=dtype)
        for i, (_, name, dtype) in enumerate(expected)
    }


def parse_charge(value: str | None, title: str) -> int:
    if value:
        digits = "".join(ch for ch in value if ch.isdigit())
        if digits:
            return int(digits)

    match = TITLE_CHARGE_RE.search(title)
    if match:
        return int(match.group(1))
    return 0


def parse_title_value(regex: re.Pattern[str], title: str, default: float = 0.0) -> float:
    match = regex.search(title)
    return float(match.group(1)) if match else default


def ensure_tof2mz_capacity(tof2mz: np.ndarray, max_tof: int) -> np.ndarray:
    if max_tof < len(tof2mz):
        return tof2mz

    new_len = len(tof2mz)
    while max_tof >= new_len:
        new_len *= 2
    grown = np.empty(new_len, dtype=np.float32)
    grown[: len(tof2mz)] = tof2mz
    grown[len(tof2mz) :] = np.nan
    return grown


def link_or_copy_pmsms(source: Path, destination: Path, copy_fragments: bool) -> None:
    if destination.exists() or destination.is_symlink():
        return
    if copy_fragments:
        shutil.copytree(source, destination, symlinks=True)
    else:
        destination.symlink_to(source.resolve(), target_is_directory=True)


def iter_mgf_spectra(path: Path):
    in_spectrum = False
    title = ""
    pepmass = None
    rt_seconds = None
    charge = None
    mz_values: list[float] = []

    with open(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line == "BEGIN IONS":
                if in_spectrum:
                    raise ValueError(f"Nested BEGIN IONS at {path}:{line_number}")
                in_spectrum = True
                title = ""
                pepmass = None
                rt_seconds = None
                charge = None
                mz_values = []
                continue
            if line == "END IONS":
                if not in_spectrum:
                    raise ValueError(f"END IONS without BEGIN IONS at {path}:{line_number}")
                if pepmass is None or rt_seconds is None:
                    raise ValueError(f"Missing PEPMASS or RTINSECONDS before {path}:{line_number}")
                yield {
                    "title": title,
                    "pepmass": pepmass,
                    "rt_seconds": rt_seconds,
                    "charge": parse_charge(charge, title),
                    "mz_values": mz_values,
                }
                in_spectrum = False
                continue

            if not in_spectrum:
                continue
            if line.startswith("TITLE="):
                title = line.removeprefix("TITLE=").strip().strip('"')
            elif line.startswith("PEPMASS="):
                pepmass = float(line.removeprefix("PEPMASS=").split()[0])
            elif line.startswith("RTINSECONDS="):
                rt_seconds = float(line.removeprefix("RTINSECONDS=").split()[0])
            elif line.startswith("CHARGE="):
                charge = line.removeprefix("CHARGE=").strip()
            elif line[0].isdigit():
                mz_values.append(float(line.split(maxsplit=1)[0]))


def convert(args: argparse.Namespace) -> dict[str, int | str]:
    source_root = args.source_root.resolve()
    mgf_path = args.mgf.resolve() if args.mgf else source_root / "plain" / "mgf.mgf"
    pmsms_dir = args.pmsms_mmappet.resolve() if args.pmsms_mmappet else source_root / "pmsms.mmappet"
    output = args.output.resolve()

    if output.exists() and not args.force:
        raise FileExistsError(f"{output} already exists; pass --force to reuse/overwrite metadata")

    output.mkdir(parents=True, exist_ok=True)
    link_or_copy_pmsms(pmsms_dir, output / "pmsms.mmappet", args.copy_fragments)

    dataindex = load_dataindex(pmsms_dir / "dataindex.mmappet")
    n_indexed = len(dataindex["precursor_idx"])
    n_expected = min(n_indexed, args.limit) if args.limit else n_indexed
    frag_tof = np.memmap(pmsms_dir / "0.bin", mode="r", dtype=np.uint32)

    tof2mz = np.empty(1024 * 1024, dtype=np.float32)
    tof2mz[:] = np.nan

    precursor_idx: list[int] = []
    mz: list[float] = []
    rt: list[float] = []
    inv_ion_mobility: list[float] = []
    charges: list[int] = []
    fragment_spectrum_start: list[int] = []
    fragment_event_cnt: list[int] = []
    manifest_rows: list[tuple[int, int, int, int, int, int]] = []

    started = time.perf_counter()
    for spectrum_0based, spectrum in enumerate(iter_mgf_spectra(mgf_path)):
        if spectrum_0based >= n_expected:
            break

        idx = int(dataindex["idx"][spectrum_0based])
        count = int(dataindex["size"][spectrum_0based])
        pidx = int(dataindex["precursor_idx"][spectrum_0based])
        mz_values = spectrum["mz_values"]

        if len(mz_values) != count:
            raise ValueError(
                f"Spectrum {spectrum_0based + 1} has {len(mz_values)} MGF peaks "
                f"but dataindex says {count}"
            )

        title_pidx_match = TITLE_PRECURSOR_RE.search(spectrum["title"])
        if title_pidx_match and int(title_pidx_match.group(1)) != pidx:
            raise ValueError(
                f"Spectrum {spectrum_0based + 1} title precursor_idx="
                f"{title_pidx_match.group(1)} but dataindex says {pidx}"
            )

        tofs = np.asarray(frag_tof[idx : idx + count], dtype=np.uint32)
        if len(tofs) != count:
            raise ValueError(f"Fragment tof table ended early at spectrum {spectrum_0based + 1}")
        if count:
            tof2mz = ensure_tof2mz_capacity(tof2mz, int(tofs.max()))
            tof2mz[tofs] = np.asarray(mz_values, dtype=np.float32)

        precursor_idx.append(pidx)
        mz.append(float(spectrum["pepmass"]))
        rt.append(float(spectrum["rt_seconds"]))
        inv_ion_mobility.append(parse_title_value(TITLE_IIM_RE, spectrum["title"]))
        charges.append(int(spectrum["charge"]))
        fragment_spectrum_start.append(idx)
        fragment_event_cnt.append(count)
        manifest_rows.append((spectrum_0based + 1, spectrum_0based, pidx, int(spectrum["charge"]), idx, count))

        if args.progress and (spectrum_0based + 1) % args.progress == 0:
            elapsed = time.perf_counter() - started
            print(
                f"[convert] {spectrum_0based + 1}/{n_expected} spectra "
                f"({elapsed:.1f}s)",
                file=sys.stderr,
            )

    if len(precursor_idx) != n_expected:
        raise ValueError(f"MGF ended after {len(precursor_idx)} spectra; expected {n_expected}")

    last_used_tof = int(np.flatnonzero(~np.isnan(tof2mz))[-1]) if np.any(~np.isnan(tof2mz)) else 0
    tof2mz = tof2mz[: last_used_tof + 1]

    tof2mz_dir = output / "tof2mz.mmappet"
    tof2mz_dir.mkdir(exist_ok=True)
    tof2mz.astype(np.float32).tofile(tof2mz_dir / "0.bin")
    (tof2mz_dir / "schema.txt").write_text("float32 x\n")
    (tof2mz_dir / "shape.txt").write_text(f"{len(tof2mz)}\n")

    df = pd.DataFrame(
        {
            "precursor_idx": np.asarray(precursor_idx, dtype=np.uint64),
            "mz": np.asarray(mz, dtype=np.float64),
            "rt": np.asarray(rt, dtype=np.float64),
            "inv_ion_mobility": np.asarray(inv_ion_mobility, dtype=np.float64),
            "charges": np.asarray(charges, dtype=np.int64),
            "fragment_spectrum_start": np.asarray(fragment_spectrum_start, dtype=np.uint64),
            "fragment_event_cnt": np.asarray(fragment_event_cnt, dtype=np.uint64),
        }
    )
    df.to_parquet(output / "precursors.parquet", index=False)

    with open(output / "manifest.tsv", "w") as handle:
        handle.write(
            "spectrum_1based\tsource_parquet_row_0based\tprecursor_idx\tcharge\t"
            "fragment_spectrum_start\tfragment_event_cnt\n"
        )
        for row in manifest_rows:
            handle.write("\t".join(str(value) for value in row) + "\n")

    return {
        "output": str(output),
        "spectra": len(precursor_idx),
        "tof2mz_len": len(tof2mz),
        "pmsms_link": str((output / "pmsms.mmappet").resolve()),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "source_root",
        type=Path,
        help="Directory containing pmsms.mmappet/ and plain/mgf.mgf.",
    )
    parser.add_argument("output", type=Path, help="Output .pmsms directory.")
    parser.add_argument("--mgf", type=Path, help="Override MGF path.")
    parser.add_argument("--pmsms-mmappet", type=Path, help="Override pmsms.mmappet path.")
    parser.add_argument("--limit", type=int, help="Convert only the first N spectra.")
    parser.add_argument("--force", action="store_true", help="Reuse an existing output directory.")
    parser.add_argument(
        "--copy-fragments",
        action="store_true",
        help="Copy pmsms.mmappet instead of symlinking it.",
    )
    parser.add_argument(
        "--progress",
        type=int,
        default=50_000,
        help="Print progress every N spectra; use 0 to disable.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.progress == 0:
        args.progress = None
    summary = convert(args)
    for key, value in summary.items():
        print(f"{key}: {value}")


if __name__ == "__main__":
    main()
