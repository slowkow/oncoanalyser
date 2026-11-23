#!/usr/bin/env python3
"""Resources used by each step of one or more oncoanalyser runs, from their Nextflow trace files.

    tools/trace_summary.py [--tsv] [--warn 0.8] <execution_trace_*.txt> ...

One row per process (the last part of the Nextflow process name), over the tasks that completed. A task
that appears in several traces (a resumed run) is counted once.

    tasks       completed tasks; "+N" adds tasks that failed or were retried
    attempt     highest attempt among completed tasks (above 1: the task was retried, often for memory)
    cpus        CPUs requested (largest)
    mem_req     memory requested (largest), GB
    peak_max    largest peak resident memory of a task, GB
    peak_p50    median peak resident memory, GB
    use         peak_max / mem_req of the same task, largest over the tasks
    wall_max    longest task, hours
    cpu_h       CPU hours: run time times CPU use, summed over tasks
    cpu_use     CPU hours over (CPUs requested times run time): how well the requested CPUs were used

Listed after the table: tasks whose peak memory reached the --warn fraction of their request, and tasks
that failed (each retry is a row; exit 137 is usually out of memory).
Columns the trace does not have (the default Nextflow trace lacks cpus, memory and attempt) are left
empty. Works with the formatted trace (`1.5 GB`, `1h 2m`) and with `trace.raw = true`.
"""
import argparse
import re
import statistics
import sys

UNITS = {"B": 1, "KB": 1 << 10, "MB": 1 << 20, "GB": 1 << 30, "TB": 1 << 40, "PB": 1 << 50}
TIME_UNITS = {"ms": 0.001, "s": 1, "m": 60, "h": 3600, "d": 86400}
GB = 1 << 30


def parse_bytes(text):
    """bytes from '1.5 GB' or a raw number; None if absent"""
    text = (text or "").strip()
    if text in ("", "-"):
        return None
    if re.fullmatch(r"[\d.]+", text):
        return float(text)
    m = re.fullmatch(r"([\d.]+)\s*([KMGTP]?B)", text)
    return float(m.group(1)) * UNITS[m.group(2)] if m else None


def parse_seconds(text, raw_unit=0.001):
    """seconds from '1h 2m 3s' or a raw number of milliseconds; None if absent"""
    text = (text or "").strip()
    if text in ("", "-"):
        return None
    if re.fullmatch(r"[\d.]+", text):
        return float(text) * raw_unit
    parts = re.findall(r"([\d.]+)(ms|s|m|h|d)", text)
    return sum(float(v) * TIME_UNITS[u] for v, u in parts) if parts else None


def parse_number(text):
    text = (text or "").strip().rstrip("%")
    try:
        return float(text)
    except ValueError:
        return None


def task_tag(row):
    """the tag (sample) of a task, from the name when the trace has no tag column"""
    if row.get("tag") not in (None, "", "-"):
        return row["tag"]
    m = re.search(r"\((.*)\)\s*$", row.get("name", ""))
    return m.group(1) if m else ""


def process_name(row):
    name = row.get("process") or row.get("name", "").split(" (")[0]
    return name.split(":")[-1]


def read_traces(paths):
    """completed tasks, one per hash; the number of other (failed, aborted) task records per process; and
    those other records"""
    done, other, failed = {}, {}, []
    for path in paths:
        with open(path) as f:
            header = f.readline().rstrip("\n").split("\t")
            for line in f:
                row = dict(zip(header, line.rstrip("\n").split("\t")))
                status = row.get("status", "")
                if status == "COMPLETED":
                    done.setdefault(row.get("hash") or row.get("task_id"), row)
                elif status not in ("CACHED", ""):
                    name = process_name(row)
                    other[name] = other.get(name, 0) + 1
                    if status == "FAILED":
                        failed.append(row)
    return list(done.values()), other, failed


def gb(value):
    return "" if value is None else f"{value / GB:.1f}"


def summarise(tasks, other):
    by_process = {}
    for row in tasks:
        by_process.setdefault(process_name(row), []).append(row)
    rows = []
    for name, group in by_process.items():
        peaks = [p for p in (parse_bytes(r.get("peak_rss")) for r in group) if p is not None]
        reqs = [m for m in (parse_bytes(r.get("memory")) for r in group) if m is not None]
        uses = []
        for r in group:
            p, m = parse_bytes(r.get("peak_rss")), parse_bytes(r.get("memory"))
            if p is not None and m:
                uses.append(p / m)
        cpus = [c for c in (parse_number(r.get("cpus")) for r in group) if c is not None]
        attempts = [a for a in (parse_number(r.get("attempt")) for r in group) if a is not None]
        walls = [w for w in (parse_seconds(r.get("realtime")) for r in group) if w is not None]
        cpu_h = 0.0
        slot_h = 0.0
        for r in group:
            wall, use, n = parse_seconds(r.get("realtime")), parse_number(r.get("%cpu")), parse_number(r.get("cpus"))
            if wall is not None and use is not None:
                cpu_h += wall * use / 100 / 3600
                if n:
                    slot_h += wall * n / 3600
        rows.append({
            "process": name,
            "tasks": f"{len(group)}" + (f"+{other[name]}" if other.get(name) else ""),
            "attempt": f"{max(attempts):.0f}" if attempts else "",
            "cpus": f"{max(cpus):.0f}" if cpus else "",
            "mem_req": gb(max(reqs)) if reqs else "",
            "peak_max": gb(max(peaks)) if peaks else "",
            "peak_p50": gb(statistics.median(peaks)) if peaks else "",
            "use": f"{max(uses):.0%}" if uses else "",
            "wall_max": f"{max(walls) / 3600:.2f}" if walls else "",
            "cpu_h": f"{cpu_h:.2f}",
            "cpu_use": f"{cpu_h / slot_h:.0%}" if slot_h else "",
            "_cpu_h": cpu_h,
        })
    rows.sort(key=lambda r: -r["_cpu_h"])
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("traces", nargs="+")
    parser.add_argument("--tsv", action="store_true", help="tab-separated output")
    parser.add_argument("--warn", type=float, default=0.8, help="list tasks whose peak reached this fraction of their memory request")
    args = parser.parse_args()

    tasks, other, failed = read_traces(args.traces)
    rows = summarise(tasks, other)
    columns = ["process", "tasks", "attempt", "cpus", "mem_req", "peak_max", "peak_p50", "use", "wall_max", "cpu_h", "cpu_use"]
    if args.tsv:
        print("\t".join(columns))
        for r in rows:
            print("\t".join(r[c] for c in columns))
    else:
        widths = {c: max(len(c), *(len(r[c]) for r in rows)) if rows else len(c) for c in columns}
        print("  ".join(c.ljust(widths[c]) if c == "process" else c.rjust(widths[c]) for c in columns))
        for r in rows:
            print("  ".join(r[c].ljust(widths[c]) if c == "process" else r[c].rjust(widths[c]) for c in columns))
        print(f"total CPU hours {sum(r['_cpu_h'] for r in rows):.1f} over {len(tasks)} tasks")

    near = []
    for r in tasks:
        p, m = parse_bytes(r.get("peak_rss")), parse_bytes(r.get("memory"))
        if p is not None and m and p / m >= args.warn:
            near.append((p / m, process_name(r), task_tag(r), p, m))
    if near:
        out = sys.stderr if args.tsv else sys.stdout
        print(f"\ntasks at or above {args.warn:.0%} of their memory request:", file=out)
        for use, name, tag, p, m in sorted(near, reverse=True):
            print(f"  {name} {tag}: peak {gb(p)} of {gb(m)} GB ({use:.0%})", file=out)
    if failed:
        out = sys.stderr if args.tsv else sys.stdout
        # exit 137: killed, usually for memory; 140: over the time limit (Slurm, some Batch setups)
        print("\ntasks that failed (each retry is a row; exit 137 is usually out of memory; aborted tasks are only counted):", file=out)
        for r in failed:
            p, m = parse_bytes(r.get("peak_rss")), parse_bytes(r.get("memory"))
            memory = f", peak {gb(p)} of {gb(m)} GB" if p is not None and m else (f", peak {gb(p)} GB" if p is not None else "")
            print(f"  {process_name(r)} {task_tag(r)}: {r.get('status', '')} exit {r.get('exit', '')}, attempt {r.get('attempt', '') or '?'}{memory}", file=out)


if __name__ == "__main__":
    main()
