"""Render an HTML telemetry report from sage + sagepy runs.

Inputs (telemetry JSONs from bench/run_with_telemetry.py — ideally invoked with
--keep-samples so time-series charts can be drawn):

  --sage    path to sage telemetry JSON
  --sagepy  path to sagepy telemetry JSON
  --output  destination HTML path

Optional PSM parity (renders identification overlap and hyperscore stats):
  --sage-psms    sage `results.sage.parquet`
  --sagepy-psms  sagepy_psms.parquet
  --q            FDR threshold (default 0.01)

Output is a single self-contained HTML file: matplotlib renders to inline
base64 PNGs, so the report has no runtime dependencies.
"""
from __future__ import annotations

import argparse
import base64
import datetime as dt
import html
import io
import json
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

SAGE_COLOR = "#1f77b4"
SAGEPY_COLOR = "#ff7f0e"
PHASE_PALETTE = [
    "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6",
    "#ffff99", "#1f78b4", "#33a02c", "#e31a1c", "#ff7f00",
]


# ---------------- IO + helpers ----------------

def load(p: Path) -> dict[str, Any]:
    with open(p) as f:
        d = json.load(f)
    # Defensive: collapse adjacent same-name phases in case the JSON was produced
    # by an older wrapper that didn't already do this.
    phases = d.get("phases", [])
    merged: list[dict[str, Any]] = []
    for ph in phases:
        if merged and merged[-1]["name"] == ph["name"]:
            prev = merged[-1]
            prev["end_t"] = ph["end_t"]
            prev["wall_s"] = round(prev["end_t"] - prev["start_t"], 3)
            for k in ("cpu_user_s", "cpu_sys_s", "minflt", "majflt", "vcsw", "ivcsw",
                      "read_bytes", "write_bytes", "samples"):
                if k in prev and k in ph:
                    prev[k] = (prev.get(k) or 0) + (ph.get(k) or 0)
            for k in ("peak_rss_kb", "peak_anon_kb", "peak_vm_size_kb",
                      "threads_p99", "threads_max"):
                if k in prev and k in ph:
                    prev[k] = max(prev.get(k) or 0, ph.get(k) or 0)
            for k in ("mean_rss_kb",):
                if k in prev and k in ph:
                    # Weighted average by sample counts (best effort).
                    a, b = prev.get("samples") or 1, ph.get("samples") or 1
                    prev[k] = round((prev[k] * a + ph[k] * b) / (a + b))
        else:
            merged.append(dict(ph))
    d["phases"] = merged
    return d


def fig_to_b64(fig: plt.Figure) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def img(b64: str, alt: str = "") -> str:
    return f'<img alt="{html.escape(alt)}" src="data:image/png;base64,{b64}" />'


def gib(kb: float | None) -> float | None:
    return None if kb is None else kb / 1024 / 1024


def safe(x: Any, fmt: str = ".2f") -> str:
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "—"
    if isinstance(x, (int, np.integer)):
        return f"{int(x):,}"
    if isinstance(x, float):
        return format(x, fmt)
    return str(x)


def ratio(num: float | None, den: float | None) -> str:
    if num is None or den in (None, 0):
        return "—"
    return f"{num / den:.2f}×"


# ---------------- tables ----------------

def headline_table(sage: dict, sp: dict) -> str:
    sage_o = sage.get("overall", {})
    sp_o = sp.get("overall", {})
    sage_cpu = sage["host"]["cpu_count"] or 1
    sp_cpu = sp["host"]["cpu_count"] or 1

    rows = [
        ("wall (s)", sage["wall_seconds"], sp["wall_seconds"], ".2f"),
        ("peak RSS (GiB)", gib(sage_o.get("peak_rss_kb")), gib(sp_o.get("peak_rss_kb")), ".2f"),
        ("mean RSS (GiB)", gib(sage_o.get("mean_rss_kb")), gib(sp_o.get("mean_rss_kb")), ".2f"),
        ("CPU user (s)", sage_o.get("cpu_user_s"), sp_o.get("cpu_user_s"), ".1f"),
        ("CPU sys (s)", sage_o.get("cpu_sys_s"), sp_o.get("cpu_sys_s"), ".1f"),
        ("parallel eff",
         (sage_o.get("cpu_user_s") or 0) / max(sage["wall_seconds"], 1e-9) / sage_cpu,
         (sp_o.get("cpu_user_s") or 0) / max(sp["wall_seconds"], 1e-9) / sp_cpu, ".2f"),
        ("threads p99", sage_o.get("threads_p99"), sp_o.get("threads_p99"), ".0f"),
        ("read bytes (GiB)",
         (sage_o.get("read_bytes") or 0) / 1024**3,
         (sp_o.get("read_bytes") or 0) / 1024**3, ".2f"),
        ("write bytes (GiB)",
         (sage_o.get("write_bytes") or 0) / 1024**3,
         (sp_o.get("write_bytes") or 0) / 1024**3, ".2f"),
        ("major faults", sage_o.get("majflt"), sp_o.get("majflt"), ",d"),
        ("minor faults", sage_o.get("minflt"), sp_o.get("minflt"), ",d"),
        ("ctx switch (vol)", sage_o.get("vcsw"), sp_o.get("vcsw"), ",d"),
        ("ctx switch (invol)", sage_o.get("ivcsw"), sp_o.get("ivcsw"), ",d"),
    ]
    parts = ["<table class='t headline'><thead><tr>",
             "<th>metric</th><th>sage</th><th>sagepy</th>",
             "<th>Δ (sagepy − sage)</th><th>ratio</th></tr></thead><tbody>"]
    for name, sv, pv, fmt in rows:
        if sv is not None and pv is not None and isinstance(sv, (int, float)) and isinstance(pv, (int, float)):
            delta = pv - sv
            r = ratio(pv, sv)
        else:
            delta = None
            r = "—"
        parts.append(
            f"<tr><td>{html.escape(name)}</td>"
            f"<td>{safe(sv, fmt)}</td>"
            f"<td>{safe(pv, fmt)}</td>"
            f"<td>{safe(delta, fmt)}</td>"
            f"<td>{html.escape(r)}</td></tr>"
        )
    parts.append("</tbody></table>")
    return "\n".join(parts)


def phase_table(sage: dict, sp: dict) -> str:
    sage_p = {p["name"]: p for p in sage["phases"]}
    sp_p = {p["name"]: p for p in sp["phases"]}
    all_names = list(dict.fromkeys(list(sage_p) + list(sp_p)))

    cols = [
        ("wall (s)", "wall_s", ".2f"),
        ("cpu_user (s)", "cpu_user_s", ".1f"),
        ("cpu_sys (s)", "cpu_sys_s", ".2f"),
        ("peak RSS (GiB)", "peak_rss_kb", ".2f"),
        ("mean RSS (GiB)", "mean_rss_kb", ".2f"),
        ("majflt", "majflt", ",d"),
        ("minflt", "minflt", ",d"),
        ("threads p99", "threads_p99", ".0f"),
    ]
    head = "<tr><th>phase</th>"
    for label, _, _ in cols:
        head += f"<th>sage {html.escape(label)}</th><th>sagepy {html.escape(label)}</th>"
    head += "</tr>"

    body = []
    for name in all_names:
        s = sage_p.get(name, {})
        p = sp_p.get(name, {})
        body.append(f"<tr><td><b>{html.escape(name)}</b></td>")
        for label, key, fmt in cols:
            sv = s.get(key)
            pv = p.get(key)
            if "RSS" in label:
                sv = gib(sv) if sv is not None else None
                pv = gib(pv) if pv is not None else None
            body.append(f"<td>{safe(sv, fmt)}</td><td>{safe(pv, fmt)}</td>")
        body.append("</tr>")

    return ("<table class='t phases'><thead>" + head + "</thead><tbody>"
            + "".join(body) + "</tbody></table>")


# ---------------- plots ----------------

def _phase_bands(ax: plt.Axes, phases: list[dict], y_lo: float, y_hi: float, alpha: float = 0.18) -> None:
    for i, ph in enumerate(phases):
        if ph["name"] == "tail":
            continue
        color = PHASE_PALETTE[i % len(PHASE_PALETTE)]
        ax.axvspan(ph["start_t"], ph["end_t"], alpha=alpha, color=color, lw=0)
        mid = (ph["start_t"] + ph["end_t"]) / 2
        ax.text(mid, y_hi, ph["name"], ha="center", va="bottom", fontsize=8,
                color="#333", clip_on=False)


def plot_phase_stack(sage: dict, sp: dict) -> str:
    """Side-by-side stacked bar of phase wall times."""
    sage_phases = [p for p in sage["phases"] if p["name"] != "tail"]
    sp_phases = [p for p in sp["phases"] if p["name"] != "tail"]
    all_names = list(dict.fromkeys([p["name"] for p in sage_phases]
                                   + [p["name"] for p in sp_phases]))

    fig, ax = plt.subplots(figsize=(10.5, 4.2))
    width = 0.55

    def bar_for(tool_phases, x):
        bottom = 0.0
        for i, name in enumerate(all_names):
            ph = next((q for q in tool_phases if q["name"] == name), None)
            if ph is None:
                continue
            color = PHASE_PALETTE[i % len(PHASE_PALETTE)]
            h = ph["wall_s"]
            ax.bar(x, h, width=width, bottom=bottom, color=color, edgecolor="white", lw=0.5)
            bottom += h

    bar_for(sage_phases, 0)
    bar_for(sp_phases, 1)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["sage", "sagepy"])
    ax.set_ylabel("wall seconds")
    ax.set_title("Phase wall-time breakdown")
    ax.grid(axis="y", alpha=0.3)
    ax.set_ylim(0, max(sage["wall_seconds"], sp["wall_seconds"]) * 1.05)
    handles = [
        Patch(facecolor=PHASE_PALETTE[i % len(PHASE_PALETTE)], edgecolor="white", label=name)
        for i, name in enumerate(all_names)
    ]
    ax.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        fontsize=8,
        title="phase",
        title_fontsize=9,
    )
    return fig_to_b64(fig)


def _series(samples: list[dict], key: str) -> tuple[np.ndarray, np.ndarray]:
    if not samples:
        return np.array([]), np.array([])
    t = np.array([s["t"] for s in samples], dtype=float)
    v = np.array([s.get(key, 0) for s in samples], dtype=float)
    return t, v


def plot_rss(sage: dict, sp: dict) -> str | None:
    s_samples = sage.get("samples")
    p_samples = sp.get("samples")
    if not s_samples or not p_samples:
        return None
    fig, axes = plt.subplots(2, 1, figsize=(10, 5.2), sharex=True)
    for ax, telem, color, samples in zip(
        axes, [sage, sp], [SAGE_COLOR, SAGEPY_COLOR], [s_samples, p_samples]
    ):
        t, rss = _series(samples, "vm_rss_kb")
        anon = _series(samples, "rss_anon_kb")[1]
        ax.plot(t, rss / 1024 / 1024, color=color, lw=1.4, label="VmRSS")
        ax.plot(t, anon / 1024 / 1024, color=color, lw=1.0, ls=":", label="RssAnon")
        _phase_bands(ax, telem["phases"], 0, max(rss / 1024 / 1024) if len(rss) else 1.0)
        ax.set_ylabel(f"{telem['tool']} GiB")
        ax.legend(loc="upper right", fontsize=8)
        ax.grid(alpha=0.3)
    axes[-1].set_xlabel("seconds since launch")
    fig.suptitle("Resident memory over time")
    return fig_to_b64(fig)


def plot_cpu(sage: dict, sp: dict) -> str | None:
    s_samples = sage.get("samples")
    p_samples = sp.get("samples")
    if not s_samples or not p_samples:
        return None

    fig, axes = plt.subplots(2, 1, figsize=(10, 5.2), sharex=True)
    for ax, telem, color, samples in zip(
        axes, [sage, sp], [SAGE_COLOR, SAGEPY_COLOR], [s_samples, p_samples]
    ):
        t, ut = _series(samples, "utime_s")
        st = _series(samples, "stime_s")[1]
        if t.size < 2:
            continue
        cum = ut + st
        # smoothed effective cores: rolling 1.0s window
        win = max(2, int(round(1.0 / max(telem["sample_interval_s"], 1e-3))))
        cores = np.zeros_like(t)
        for i in range(t.size):
            j = max(0, i - win)
            dt = t[i] - t[j]
            cores[i] = (cum[i] - cum[j]) / dt if dt > 0 else 0.0
        ax.plot(t, cores, color=color, lw=1.2)
        ax.axhline(telem["host"]["cpu_count"], color="#888", lw=0.8, ls="--",
                   label=f"{telem['host']['cpu_count']} cores")
        _phase_bands(ax, telem["phases"], 0, max(cores) if cores.size else 1.0)
        ax.set_ylabel(f"{telem['tool']} effective cores")
        ax.legend(loc="upper right", fontsize=8)
        ax.grid(alpha=0.3)
    axes[-1].set_xlabel("seconds since launch")
    fig.suptitle("CPU utilization (1s rolling, effective cores)")
    return fig_to_b64(fig)


def plot_io_faults(sage: dict, sp: dict) -> str | None:
    s_samples = sage.get("samples")
    p_samples = sp.get("samples")
    if not s_samples or not p_samples:
        return None
    fig, axes = plt.subplots(2, 2, figsize=(11, 5.5), sharex="col")

    def rate(t: np.ndarray, c: np.ndarray, win: int) -> np.ndarray:
        if t.size < 2:
            return np.zeros_like(t)
        out = np.zeros_like(t)
        for i in range(t.size):
            j = max(0, i - win)
            dt = t[i] - t[j]
            out[i] = (c[i] - c[j]) / dt if dt > 0 else 0.0
        return out

    for col_idx, (telem, color, samples) in enumerate(
        [(sage, SAGE_COLOR, s_samples), (sp, SAGEPY_COLOR, p_samples)]
    ):
        win = max(2, int(round(1.0 / max(telem["sample_interval_s"], 1e-3))))
        t, mn = _series(samples, "minflt")
        mj = _series(samples, "majflt")[1]
        rb = _series(samples, "read_bytes")[1]
        wb = _series(samples, "write_bytes")[1]

        ax = axes[0, col_idx]
        ax.plot(t, rate(t, mn, win) / 1000, color=color, lw=1.1, label="minor /1k")
        ax.plot(t, rate(t, mj, win),       color=color, lw=1.1, ls="--", label="major")
        _phase_bands(ax, telem["phases"], 0, max(rate(t, mn, win) / 1000) if t.size else 1.0)
        ax.set_ylabel(f"{telem['tool']} faults/s")
        ax.legend(loc="upper right", fontsize=8)
        ax.grid(alpha=0.3)

        ax = axes[1, col_idx]
        ax.plot(t, rate(t, rb, win) / 1024 / 1024, color=color, lw=1.1, label="read MiB/s")
        ax.plot(t, rate(t, wb, win) / 1024 / 1024, color=color, lw=1.1, ls="--", label="write MiB/s")
        _phase_bands(ax, telem["phases"], 0, max(rate(t, rb, win) / 1024 / 1024) if t.size else 1.0)
        ax.set_ylabel(f"{telem['tool']} MiB/s")
        ax.set_xlabel("seconds since launch")
        ax.legend(loc="upper right", fontsize=8)
        ax.grid(alpha=0.3)

    fig.suptitle("Page faults and disk IO (1s rolling rate)")
    return fig_to_b64(fig)


# ---------------- PSM parity (lightweight, derived from compare_psms.py) ----------------

def parity_section(sage_psms: Path | None, sagepy_psms: Path | None, q: float) -> str:
    if not sage_psms or not sagepy_psms:
        return ""
    sp = pd.read_parquet(sagepy_psms).rename(columns={
        "decoy": "is_decoy",
        "sequence_modified": "peptide",
        "sequence": "stripped_peptide",
    })
    sg = pd.read_parquet(sage_psms)

    def topf(df):
        return df[(df["rank"] == 1) & (df["is_decoy"] == False) & (df["spectrum_q"] <= q)]

    sp_top = topf(sp)
    sg_top = topf(sg)
    sp_peps = set(sp_top["peptide"])
    sg_peps = set(sg_top["peptide"])
    sp_strip = set(sp_top["stripped_peptide"])
    sg_strip = set(sg_top["stripped_peptide"])

    rows = [
        ("raw PSMs",            len(sp), len(sg)),
        ("rank-1",              int((sp["rank"] == 1).sum()), int((sg["rank"] == 1).sum())),
        ("targets",             int((~sp["is_decoy"].astype(bool)).sum()),
                                int((~sg["is_decoy"].astype(bool)).sum())),
        (f"top1 target q<{q}",  len(sp_top), len(sg_top)),
        ("unique peptides (mod)",     len(sp_peps), len(sg_peps)),
        ("shared peptides (mod)",     len(sp_peps & sg_peps), "—"),
        ("sagepy-only peptides (mod)", len(sp_peps - sg_peps), "—"),
        ("sage-only peptides (mod)",   "—", len(sg_peps - sp_peps)),
        ("unique peptides (stripped)",  len(sp_strip), len(sg_strip)),
        ("shared peptides (stripped)",  len(sp_strip & sg_strip), "—"),
        ("hyperscore mean (top1 q<1%)",
         f"{sp_top['hyperscore'].mean():.2f}",
         f"{sg_top['hyperscore'].mean():.2f}"),
        ("hyperscore median (top1 q<1%)",
         f"{sp_top['hyperscore'].median():.2f}",
         f"{sg_top['hyperscore'].median():.2f}"),
    ]
    parts = ["<h2>PSM parity</h2><table class='t parity'><thead>"
             "<tr><th>metric</th><th>sagepy</th><th>sage</th></tr></thead><tbody>"]
    for name, sv, gv in rows:
        parts.append(
            f"<tr><td>{html.escape(name)}</td>"
            f"<td>{safe(sv, ',d') if isinstance(sv, int) else html.escape(str(sv))}</td>"
            f"<td>{safe(gv, ',d') if isinstance(gv, int) else html.escape(str(gv))}</td></tr>"
        )
    parts.append("</tbody></table>")
    return "\n".join(parts)


# ---------------- top-level render ----------------

CSS = """
body { font: 14px/1.45 -apple-system, Segoe UI, Roboto, sans-serif; max-width: 1100px;
       margin: 24px auto; padding: 0 16px; color: #1a1a1a; }
h1 { font-size: 22px; margin-bottom: 4px; }
h2 { font-size: 17px; margin-top: 28px; border-bottom: 1px solid #ddd; padding-bottom: 4px; }
h3 { font-size: 14px; margin-top: 16px; color: #555; }
.meta { color: #555; font-size: 12px; }
.meta dl { display: grid; grid-template-columns: max-content 1fr; gap: 2px 12px; margin: 8px 0; }
.meta dt { font-weight: 600; }
.meta dd { margin: 0; font-family: ui-monospace, Menlo, Consolas, monospace; }
table.t { border-collapse: collapse; font-size: 12.5px; margin: 8px 0; }
table.t th, table.t td { padding: 4px 10px; border-bottom: 1px solid #eee; text-align: right; }
table.t th { background: #f6f6f6; text-align: right; }
table.t td:first-child, table.t th:first-child { text-align: left; }
table.t tbody tr:hover { background: #fafafa; }
img { max-width: 100%; height: auto; display: block; margin: 4px 0 16px; }
.figcap { color: #777; font-size: 11.5px; margin-top: -10px; margin-bottom: 16px; }
"""


def metadata_section(sage: dict, sp: dict) -> str:
    def cmd_line(d: dict) -> str:
        return " ".join(html.escape(c) for c in d.get("command", []))
    when = dt.datetime.now().isoformat(timespec="seconds")
    return f"""
<section class="meta">
<dl>
  <dt>generated</dt><dd>{when}</dd>
  <dt>host</dt><dd>{html.escape(sage["host"]["uname"])}</dd>
  <dt>cpu_count</dt><dd>{sage["host"]["cpu_count"]}</dd>
  <dt>sage label</dt><dd>{html.escape(sage.get("label") or "—")}</dd>
  <dt>sagepy label</dt><dd>{html.escape(sp.get("label") or "—")}</dd>
  <dt>sage cmd</dt><dd>{cmd_line(sage)}</dd>
  <dt>sagepy cmd</dt><dd>{cmd_line(sp)}</dd>
  <dt>sage exit</dt><dd>{sage.get("exit_code")}</dd>
  <dt>sagepy exit</dt><dd>{sp.get("exit_code")}</dd>
</dl>
</section>
"""


def render(sage: dict, sp: dict, sage_psms: Path | None,
           sagepy_psms: Path | None, q: float) -> str:
    rss_b64 = plot_rss(sage, sp)
    cpu_b64 = plot_cpu(sage, sp)
    iof_b64 = plot_io_faults(sage, sp)
    phase_b64 = plot_phase_stack(sage, sp)

    rss_html = (img(rss_b64, "rss") + "<p class='figcap'>"
                "Solid = VmRSS (resident set), dotted = RssAnon (private anonymous). "
                "Phase shading per tool.</p>") if rss_b64 else \
        "<p class='figcap'>RSS time-series unavailable (run with <code>--keep-samples</code>).</p>"
    cpu_html = (img(cpu_b64, "cpu") + "<p class='figcap'>"
                "Effective cores in use, computed as Δ(utime+stime)/Δt over a 1 s window. "
                "Dashed line = total cores on host.</p>") if cpu_b64 else ""
    iof_html = img(iof_b64, "faults+io") if iof_b64 else ""

    parity = parity_section(sage_psms, sagepy_psms, q)

    return f"""<!doctype html>
<html><head><meta charset="utf-8"><title>sage vs sagepy telemetry</title>
<style>{CSS}</style></head><body>
<h1>sage vs sagepy — telemetry report</h1>
{metadata_section(sage, sp)}

<h2>Headline</h2>
{headline_table(sage, sp)}

<h2>Phase wall-time breakdown</h2>
{img(phase_b64, "phase stack")}

<h2>Per-phase metrics</h2>
{phase_table(sage, sp)}

<h2>Memory over time</h2>
{rss_html}

<h2>CPU utilization over time</h2>
{cpu_html}

<h2>Page faults and disk IO</h2>
{iof_html}

{parity}

</body></html>
"""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sage", type=Path, required=True)
    ap.add_argument("--sagepy", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--sage-psms", type=Path, default=None)
    ap.add_argument("--sagepy-psms", type=Path, default=None)
    ap.add_argument("--q", type=float, default=0.01)
    args = ap.parse_args()

    sage = load(args.sage)
    sp = load(args.sagepy)
    html_doc = render(sage, sp, args.sage_psms, args.sagepy_psms, args.q)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(html_doc)
    print(f"[render] {args.output}  ({len(html_doc):,} bytes)")


if __name__ == "__main__":
    main()
