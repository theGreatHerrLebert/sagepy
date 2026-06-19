"""External telemetry wrapper for the sage / sagepy benchmarks.

Spawns a subprocess, polls /proc/<pid>/{status,stat,io} every 50 ms, captures
stderr/stdout line-by-line with timestamps, then partitions the sample stream
into phases by parsing each tool's existing log markers (no source patches).

Per-phase metrics emitted: wall, cpu_user, cpu_sys, peak/mean RSS, peak anon,
voluntary/involuntary context switches, minor/major page faults, IO bytes,
peak thread count.

Usage:
    python bench/run_with_telemetry.py --tool sage --output sage_telemetry.json -- \\
        /path/to/sage --fasta human.fasta -o results/sage --parquet \\
        bench/sage_config.json input.pmsms

    python bench/run_with_telemetry.py --tool sagepy --output sagepy_telemetry.json -- \\
        bench/.venv/bin/python bench/run_sagepy.py --rescore --input input.pmsms
"""
from __future__ import annotations

import argparse
import atexit
import json
import os
import re
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import Any, Callable

SAMPLE_INTERVAL = 0.05
CLK_TCK = os.sysconf(os.sysconf_names["SC_CLK_TCK"])
PAGE_SIZE = os.sysconf(os.sysconf_names["SC_PAGESIZE"])


def _read_status(pid: int) -> dict[str, int]:
    out: dict[str, int] = {}
    try:
        with open(f"/proc/{pid}/status") as f:
            for line in f:
                if line.startswith((
                    "VmRSS:", "VmHWM:", "VmPeak:", "RssAnon:", "VmSize:",
                    "voluntary_ctxt_switches:", "nonvoluntary_ctxt_switches:",
                    "Threads:",
                )):
                    name, _, rest = line.partition(":")
                    parts = rest.strip().split()
                    if parts:
                        try:
                            out[name] = int(parts[0])
                        except ValueError:
                            pass
    except FileNotFoundError:
        pass
    return out


def _read_stat(pid: int) -> dict[str, int] | None:
    try:
        with open(f"/proc/{pid}/stat") as f:
            data = f.read()
    except FileNotFoundError:
        return None
    rparen = data.rfind(")")
    fields = data[rparen + 2:].split()
    return {
        "minflt": int(fields[7]),
        "majflt": int(fields[9]),
        "utime_ticks": int(fields[11]),
        "stime_ticks": int(fields[12]),
        "num_threads": int(fields[17]),
    }


def _read_io(pid: int) -> dict[str, int]:
    out: dict[str, int] = {}
    try:
        with open(f"/proc/{pid}/io") as f:
            for line in f:
                k, _, v = line.partition(":")
                try:
                    out[k.strip()] = int(v.strip())
                except ValueError:
                    pass
    except (FileNotFoundError, PermissionError):
        pass
    return out


def sample(pid: int) -> dict[str, Any] | None:
    status = _read_status(pid)
    stat = _read_stat(pid)
    if not status or stat is None:
        return None
    io = _read_io(pid)
    return {
        "vm_rss_kb": status.get("VmRSS", 0),
        "vm_hwm_kb": status.get("VmHWM", 0),
        "vm_peak_kb": status.get("VmPeak", 0),
        "rss_anon_kb": status.get("RssAnon", 0),
        "vm_size_kb": status.get("VmSize", 0),
        "threads": status.get("Threads", stat["num_threads"]),
        "vcsw": status.get("voluntary_ctxt_switches", 0),
        "ivcsw": status.get("nonvoluntary_ctxt_switches", 0),
        "minflt": stat["minflt"],
        "majflt": stat["majflt"],
        "utime_s": stat["utime_ticks"] / CLK_TCK,
        "stime_s": stat["stime_ticks"] / CLK_TCK,
        "read_bytes": io.get("read_bytes", 0),
        "write_bytes": io.get("write_bytes", 0),
    }


SAGEPY_PHASE_RE = re.compile(r"^\[(read|digest|process|score|rescore|lda|fdr|write)\]\s")


def parse_sagepy(line: str) -> str | None:
    m = SAGEPY_PHASE_RE.match(line)
    if not m:
        return None
    phase = m.group(1)
    # Phase-end lines always end with "<float>s" (some say "in Xs", fdr just "Xs").
    if re.search(r"[0-9.]+s\s*$", line):
        return phase
    return None


SAGE_PHASE_END_PATTERNS: list[tuple[str, re.Pattern[str]]] = [
    ("digest",  re.compile(r"generated \d+ fragments, \d+ peptides in")),
    ("read",    re.compile(r"- file IO:\s+\d+ ms")),
    ("score",   re.compile(r"- search:\s+\d+ ms")),
    ("rescore", re.compile(r"fit mobility model|fit retention time model")),
    ("fdr",     re.compile(r"discovered \d+ target peptide-spectrum matches")),
    ("write",   re.compile(r"finished in \d+s")),
]


def parse_sage(line: str) -> str | None:
    for name, pat in SAGE_PHASE_END_PATTERNS:
        if pat.search(line):
            return name
    return None


PARSERS: dict[str, Callable[[str], str | None]] = {
    "sage": parse_sage,
    "sagepy": parse_sagepy,
}


def aggregate(samples: list[tuple[float, dict[str, Any]]],
              start_t: float, end_t: float) -> dict[str, Any]:
    in_phase = [s for t, s in samples if start_t <= t <= end_t]
    if not in_phase:
        return {"samples": 0, "wall_s": round(end_t - start_t, 3)}
    rss = [s["vm_rss_kb"] for s in in_phase]
    threads = sorted(s["threads"] for s in in_phase)
    first, last = in_phase[0], in_phase[-1]
    p99_idx = max(0, int(len(threads) * 0.99) - 1)
    return {
        "samples": len(in_phase),
        "wall_s": round(end_t - start_t, 3),
        "cpu_user_s": round(last["utime_s"] - first["utime_s"], 3),
        "cpu_sys_s":  round(last["stime_s"] - first["stime_s"], 3),
        "peak_rss_kb": max(rss),
        "mean_rss_kb": round(sum(rss) / len(rss)),
        "peak_anon_kb": max(s["rss_anon_kb"] for s in in_phase),
        "peak_vm_size_kb": max(s["vm_size_kb"] for s in in_phase),
        "minflt": last["minflt"] - first["minflt"],
        "majflt": last["majflt"] - first["majflt"],
        "vcsw": last["vcsw"] - first["vcsw"],
        "ivcsw": last["ivcsw"] - first["ivcsw"],
        "read_bytes": last["read_bytes"] - first["read_bytes"],
        "write_bytes": last["write_bytes"] - first["write_bytes"],
        "threads_p99": threads[p99_idx],
        "threads_max": threads[-1],
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tool", choices=list(PARSERS), required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--sample-interval", type=float, default=SAMPLE_INTERVAL)
    ap.add_argument("--label", default=None)
    ap.add_argument("--keep-samples", action="store_true",
                    help="Include the full per-sample timeseries in the JSON output.")
    ap.add_argument("command", nargs=argparse.REMAINDER)
    args = ap.parse_args()

    cmd = args.command
    if cmd and cmd[0] == "--":
        cmd = cmd[1:]
    if not cmd:
        ap.error("missing command after --")

    parser = PARSERS[args.tool]

    # Defeat fully-buffered stdout when the child is Python — otherwise phase-boundary
    # log lines would all arrive at process exit and our timestamps would be wrong.
    env = os.environ.copy()
    env.setdefault("PYTHONUNBUFFERED", "1")

    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        bufsize=1,
        text=True,
        start_new_session=True,  # own process group so we can kill the whole tree
        env=env,
    )
    pid = proc.pid
    pgid = os.getpgid(pid)
    start = time.perf_counter()

    def _kill_child(*_args: Any) -> None:
        if proc.poll() is None:
            try:
                os.killpg(pgid, signal.SIGTERM)
            except ProcessLookupError:
                pass

    atexit.register(_kill_child)
    for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
        try:
            signal.signal(sig, lambda *_: (_kill_child(), sys.exit(130)))
        except (ValueError, OSError):
            pass

    samples: list[tuple[float, dict[str, Any]]] = []
    events: list[tuple[float, str, str]] = []
    log_lines: list[tuple[float, str]] = []
    stop = threading.Event()

    def sampler() -> None:
        while not stop.is_set():
            t = time.perf_counter() - start
            s = sample(pid)
            if s is None:
                break
            samples.append((t, s))
            stop.wait(args.sample_interval)

    def reader() -> None:
        assert proc.stdout is not None
        for line in proc.stdout:
            t = time.perf_counter() - start
            line = line.rstrip("\n")
            log_lines.append((t, line))
            sys.stderr.write(line + "\n")
            sys.stderr.flush()
            phase = parser(line)
            if phase:
                events.append((t, phase, line))

    th_sampler = threading.Thread(target=sampler, daemon=True)
    th_reader = threading.Thread(target=reader, daemon=True)
    th_sampler.start()
    th_reader.start()

    rc = proc.wait()
    stop.set()
    th_reader.join(timeout=5)
    th_sampler.join(timeout=2)
    end = time.perf_counter() - start

    # Collapse adjacent same-name boundary events (e.g. sage emits both retention
    # and mobility model fits as "rescore" — we want one rescore phase, not two).
    collapsed: list[tuple[float, str, str]] = []
    for ev in events:
        if collapsed and collapsed[-1][1] == ev[1]:
            collapsed[-1] = ev
        else:
            collapsed.append(ev)

    phases: list[dict[str, Any]] = []
    prev_end = 0.0
    for t, name, _line in collapsed:
        phases.append({"name": name, "start_t": round(prev_end, 3),
                       "end_t": round(t, 3), **aggregate(samples, prev_end, t)})
        prev_end = t
    if prev_end < end:
        phases.append({"name": "tail", "start_t": round(prev_end, 3),
                       "end_t": round(end, 3), **aggregate(samples, prev_end, end)})

    overall = aggregate(samples, 0.0, end)

    summary: dict[str, Any] = {
        "tool": args.tool,
        "label": args.label,
        "command": cmd,
        "exit_code": rc,
        "wall_seconds": round(end, 3),
        "n_samples": len(samples),
        "sample_interval_s": args.sample_interval,
        "host": {
            "uname": " ".join(os.uname()),
            "cpu_count": os.cpu_count(),
            "page_size": PAGE_SIZE,
            "clk_tck": CLK_TCK,
        },
        "overall": overall,
        "phases": phases,
        "events": [{"t": round(t, 3), "phase_end": name, "line": line}
                   for t, name, line in events],
    }
    if args.keep_samples:
        summary["samples"] = [{"t": round(t, 3), **s} for t, s in samples]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\n[telemetry] wrote {args.output} "
          f"({len(samples)} samples, {len(phases)} phases, exit={rc})",
          file=sys.stderr)
    sys.exit(rc)


if __name__ == "__main__":
    main()
