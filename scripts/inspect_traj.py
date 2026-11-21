# namd-aggregation-md -- NAMD small-molecule self-assembly analysis
# Copyright (C) 2026 huyang-md
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Audit trajectory files in a case directory against their logs.

Answers, in one table: how many frames each .dcd really holds, which steps they
cover, at what frame spacing, whether the file is truncated, and whether the log
ended cleanly. Then merges the segments per system so gaps and overlaps show up.

The DCD header is authoritative: NAMD rewrites NSET/NSTEP as it goes, so a file
whose header claims N frames but holds fewer was cut short, and a file whose
header claims fewer than the log's last ENERGY step stopped before that step.

    python inspect_traj.py --dir <dir>
    python inspect_traj.py --dir <dir> --json
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import namdlib  # noqa: E402
import scaffold  # noqa: E402

ENERGY_RE = re.compile(r"^ENERGY:\s+(\d+)", re.M)


def log_status(path):
    if not os.path.exists(path):
        return None, None, "no log"
    try:
        with open(path, "r", errors="replace") as fh:
            fh.seek(max(0, os.path.getsize(path) - 400_000))
            tail = fh.read()
    except OSError:
        return None, None, "unreadable"
    steps = ENERGY_RE.findall(tail)
    last = int(steps[-1]) if steps else None
    if "End of program" in tail:
        status = "clean"
    elif "FATAL ERROR" in tail:
        fatal = tail.rsplit("FATAL ERROR", 1)[1].splitlines()[0].strip()
        status = f"FATAL: {fatal[:60]}"
    else:
        status = "no end marker"
    return last, status, None


def timestep_for(directory, default=2.0):
    case_path = os.path.join(directory, scaffold.CASE_FILE)
    if os.path.exists(case_path):
        try:
            return float(scaffold.load_case(directory)["timestep_fs"])
        except Exception:  # noqa: BLE001
            pass
    for name in ("04_prod.conf", "product.conf", "03_equil_npt.conf"):
        p = os.path.join(directory, name)
        if os.path.exists(p):
            for line in open(p, encoding="utf-8", errors="replace"):
                s = line.split("#", 1)[0].split()
                if len(s) >= 2 and s[0] == "timestep":
                    try:
                        return float(s[1])
                    except ValueError:
                        pass
    return default


def row_for(path, ts_fs):
    try:
        dcd = namdlib.Dcd(path)
    except Exception as exc:  # noqa: BLE001
        return dict(file=os.path.basename(path), error=str(exc))
    dt_ps = ts_fs / 1000.0
    ps_per_frame = dcd.nsavc * dt_ps
    ns_per_frame = ps_per_frame / 1000.0
    last_step = dcd.istart + (dcd.frames - 1) * dcd.nsavc if dcd.frames else None
    log_path = path[:-4] + ".log" if path.endswith(".dcd") else None
    log_last, status, note = (None, None, "backup copy")
    if log_path and not path.endswith(".BAK"):
        log_last, status, note = log_status(log_path)

    flags = []
    if dcd.frames != dcd.nset:
        flags.append(f"header NSET={dcd.nset} but disk holds {dcd.frames}")
    if dcd.partial_bytes:
        flags.append(f"{dcd.partial_bytes} trailing bytes (partial frame)")
    if log_last is not None and last_step is not None and last_step > log_last + 1:
        flags.append(f"trajectory reaches step {last_step:,}, past log's last ENERGY {log_last:,}")
    if log_last is not None and last_step is not None and log_last > last_step + dcd.nsavc:
        flags.append(f"log reached step {log_last:,} but trajectory stops at {last_step:,}")
    if status and status.startswith("FATAL"):
        flags.append(status)
    out = dict(
        file=os.path.basename(path),
        natoms=dcd.natoms,
        frames_hdr=dcd.nset,
        frames_disk=dcd.frames,
        first_step=dcd.istart,
        last_step=last_step,
        ps_per_frame=ps_per_frame,
        span_ns=dcd.frames * ns_per_frame if dcd.frames else 0.0,
        start_ns=(dcd.istart * dt_ps / 1000.0),
        # coverage ends at the time of the last frame actually written, so a gap
        # is measured between real data points, not between run boundaries
        end_ns=(last_step * dt_ps / 1000.0) if last_step else 0.0,
        log_last=log_last,
        log_status=status,
        size=dcd.size,
        flags="; ".join(flags),
    )
    dcd.close()
    return out


PROD_RE = re.compile(r"prod", re.I)


def merge_segments(rows, include_all=False):
    """Merge production trajectories per system into covered time spans.

    Equilibration files restart their own step counter, so mixing them into a
    production timeline is meaningless; they are skipped unless include_all is
    set. Backup copies (.BAK) are always skipped.
    """
    by_n = {}
    for r in rows:
        if "error" in r:
            continue
        if r["file"].endswith(".BAK"):
            continue
        if not include_all and not PROD_RE.search(r["file"]):
            continue
        by_n.setdefault(r["natoms"], []).append(r)
    report = []
    for natoms, group in by_n.items():
        good = group
        good.sort(key=lambda g: g["first_step"])
        segs, cursor, gaps, overlaps = [], None, [], []
        spacings = sorted({g["ps_per_frame"] for g in good})
        for g in good:
            s, e = g["start_ns"], g["end_ns"]
            if cursor is None:
                segs.append([s, e])
                cursor = e
                continue
            if s > cursor + 1e-6:
                gaps.append((cursor, s))
                segs.append([s, e])
                cursor = e
            elif s <= cursor + 1e-6:
                if e > cursor:
                    overlaps.append((s, min(e, cursor)))
                    segs[-1][1] = max(segs[-1][1], e)
                    cursor = e
        covered = sum(e - s for s, e in segs)
        report.append(dict(natoms=natoms, files=len(good), segments=segs,
                           covered_ns=covered, gaps=gaps, overlaps=overlaps,
                           spacings=spacings))
    return report


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dir", required=True)
    p.add_argument("--timestep-fs", type=float, help="override the timestep (default: from conf/case)")
    p.add_argument("--all", action="store_true",
                   help="also merge equilibration trajectories into the timeline")
    p.add_argument("--json", action="store_true")
    return p.parse_args(argv)


def main(argv=None):
    a = parse_args(argv)
    d = os.path.abspath(a.dir)
    ts = a.timestep_fs or timestep_for(d)
    files = sorted(glob.glob(os.path.join(d, "*.dcd")) + glob.glob(os.path.join(d, "*.dcd.BAK")))
    files = [f for f in files if not os.path.basename(f).startswith("#")]
    if not files:
        sys.exit(f"no .dcd files in {d}")

    rows = [row_for(f, ts) for f in files]
    merged = merge_segments(rows, a.all)

    if a.json:
        print(json.dumps(dict(directory=d, timestep_fs=ts, files=rows, segments=merged),
                         indent=2, ensure_ascii=False))
        return 0

    print(f"{d}\ntimestep {ts:g} fs\n")
    hdr = (f"{'file':<24}{'kind':<7}{'natoms':>8}{'frames':>7}{'hdr':>6}"
           f"{'first':>12}{'last':>12}{'ps/fr':>8}{'ns':>8}")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        if "error" in r:
            print(f"{r['file']:<24}  ERROR: {r['error']}")
            continue
        kind = "bak" if r["file"].endswith(".BAK") else ("prod" if PROD_RE.search(r["file"]) else "equil")
        print(f"{r['file']:<24}{kind:<7}{r['natoms']:>8}{r['frames_disk']:>7}{r['frames_hdr']:>6}"
              f"{r['first_step']:>12,}{r['last_step']:>12,}{r['ps_per_frame']:>8.1f}"
              f"{r['span_ns']:>8.1f}")
    notes = [r for r in rows if r.get("flags")]
    if notes:
        print("\nnotes:")
        for r in notes:
            print(f"  {r['file']:<22} {r['flags']}")

    scope = "all trajectories" if a.all else "production trajectories"
    print(f"\ncoverage per system ({scope}, backups skipped):")
    for m in merged:
        print(f"\n  natoms {m['natoms']:,}   {m['files']} file(s)   "
              f"{m['covered_ns']:.1f} ns covered   "
              f"spacing {', '.join(f'{s:g}' for s in m['spacings'])} ps")
        for s, e in m["segments"]:
            print(f"    {s:7.2f} -> {e:7.2f} ns")
        for g0, g1 in m["gaps"]:
            print(f"    GAP      {g0:.2f} -> {g1:.2f} ns ({g1 - g0:.2f} ns missing)")
        for o0, o1 in m["overlaps"]:
            print(f"    overlap  {o0:.2f} -> {o1:.2f} ns (same time covered twice)")
        if len(m["spacings"]) > 1:
            print("    WARN: mixed frame spacing in one system; "
                  "weight by time, not by frame, when averaging")
    return 0


if __name__ == "__main__":
    sys.exit(main())
