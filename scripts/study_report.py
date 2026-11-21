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

"""One-page status of the study: can each claim be written yet?

Combines three things that must agree before a substituent comparison means anything:
  * system comparability (concentration, frame spacing, length, force-field knobs)
  * what the analysis actually measured (J-type fraction, interplanar distance, clusters)
  * whether the assembly has stopped changing

    python study_report.py --root <parent dir>
    python study_report.py --root <parent dir> --analysis-root <analysis root>
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_systems  # noqa: E402

OK, BAD, UNK = "PASS", "BLOCKED", "UNKNOWN"

RE_SLIP = re.compile(r"slip angle: mean ([\d.]+) deg, J-type \(> \d+ deg\) ([\d.]+)%")
RE_STACK = re.compile(r"stacked pairs .*?: interplanar d median ([\d.]+) A, mean ([\d.]+) A, "
                      r"([\d.]+)% within 4 A")
RE_NEAR = re.compile(r"nearest neighbour: median ([\d.]+) A, mean ([\d.]+) A")
RE_CLUST = re.compile(r"clusters \(cutoff ([\d.]+) A, (\w+)\): ([\d.]+)% of frames have "
                      r"clusters, largest ever (\d+) molecules")
RE_TREND = re.compile(r"assembly trend: mean cluster ([\d.]+) -> ([\d.]+) molecules \((\w+)\)")
RE_SANITY = re.compile(r"SANITY FAIL")


def read_metrics(analysis_dir):
    path = os.path.join(analysis_dir, "summary.txt")
    if not os.path.exists(path):
        return None
    txt = open(path, encoding="utf-8", errors="replace").read()
    m = {}
    if (hit := RE_SLIP.search(txt)):
        m["slip_mean"] = float(hit.group(1))
        m["j_type_pct"] = float(hit.group(2))
    if (hit := RE_STACK.search(txt)):
        m["d_median"] = float(hit.group(1))
        m["d_mean"] = float(hit.group(2))
        m["d_within4_pct"] = float(hit.group(3))
    if (hit := RE_NEAR.search(txt)):
        m["nearest_median"] = float(hit.group(1))
    if (hit := RE_CLUST.search(txt)):
        m["cluster_cutoff"] = float(hit.group(1))
        m["cluster_mode"] = hit.group(2)
        m["cluster_largest"] = int(hit.group(4))
    if (hit := RE_TREND.search(txt)):
        m["cluster_first"] = float(hit.group(1))
        m["cluster_last"] = float(hit.group(2))
        m["trend"] = hit.group(3)
    m["sanity_failed"] = bool(RE_SANITY.search(txt))
    return m


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--root", required=True, help="parent dir with one system per subdir")
    p.add_argument("--analysis-root", help="dir holding one analysis/ per system "
                                           "(default: <root>/../_analysis/compare)")
    p.add_argument("--resname", default="LIG")
    p.add_argument("--json", action="store_true")
    a = p.parse_args(argv)

    root = os.path.abspath(a.root)
    ana_root = a.analysis_root or os.path.join(os.path.dirname(root), "_analysis", "compare")

    systems, metrics = [], {}
    for name in sorted(os.listdir(root)):
        d = os.path.join(root, name)
        if not os.path.isdir(d):
            continue
        try:
            info = compare_systems.analyse_case(d, a.resname)
        except Exception as exc:  # noqa: BLE001
            print(f"[skip] {name}: {exc}")
            continue
        if not info:
            continue
        systems.append(info)
        metrics[info["name"]] = read_metrics(os.path.join(ana_root, info["name"]))

    if not systems:
        sys.exit(f"no systems with a PSF under {root}")

    issues = compare_systems.compare(systems)
    failed = {k for level, k, _ in issues if level == "FAIL"}

    # ---- per-system table ---------------------------------------------------
    hdr = (f"{'system':<12}{'conc mM':>9}{'prod ns':>9}{'J-type%':>9}{'d med':>8}"
           f"{'nearest':>9}{'cluster':>9}{'trend':>9}")
    print(hdr)
    print("-" * len(hdr))
    for s in systems:
        m = metrics.get(s["name"]) or {}
        print(f"{s['name']:<12}"
              f"{s.get('conc_mM', '-'):>9}"
              f"{s.get('prod_ns', '-'):>9}"
              f"{m.get('j_type_pct', '-'):>9}"
              f"{m.get('d_median', '-'):>8}"
              f"{m.get('nearest_median', '-'):>9}"
              f"{m.get('cluster_largest', '-'):>9}"
              f"{m.get('trend', '-'):>9}")

    # ---- claims -------------------------------------------------------------
    claims = []

    if any(metrics.values()) and all(metrics.get(s["name"]) for s in systems):
        jt = {s["name"]: metrics[s["name"]].get("j_type_pct") for s in systems}
        bad = [k for k, v in jt.items() if v is not None and v < 50.0]
        claims.append(("J-type slip stacking exists", OK if not bad else BAD,
                       "all systems > 50% J-type" if not bad
                       else f"below 50% in {', '.join(bad)}"))
        dm = {s["name"]: metrics[s["name"]].get("d_median") for s in systems}
        over = [k for k, v in dm.items() if v is not None and v > 4.0]
        claims.append(("Interplanar distance <= 4 A", OK if not over else BAD,
                       "median d <= 4.0 A in every system" if not over
                       else f"median > 4.0 A in {', '.join(over)}"))
        growing = [s["name"] for s in systems
                   if (metrics[s["name"]] or {}).get("trend") == "growing"]
        claims.append(("Assembly has equilibrated", OK if not growing else BAD,
                       "cluster size is stable" if not growing
                       else f"still growing in {', '.join(growing)}; size comparison premature"))
        sanity = [s["name"] for s in systems if (metrics[s["name"]] or {}).get("sanity_failed")]
        if sanity:
            claims.append(("Analysis geometry is physical", BAD,
                           f"sanity check failed in {', '.join(sanity)}"))
    else:
        claims.append(("J-type / distance measured", UNK,
                       f"no analysis summary under {ana_root}"))

    blockers = [k for k in ("concentration", "spacing", "length", "oneFourScaling")
                if k in failed]
    if blockers:
        claims.append(("Substituent comparison", BAD,
                       "systems differ in " + ", ".join(blockers)
                       + " - fix those before ranking aggregation"))
    else:
        claims.append(("Substituent comparison", OK,
                       "concentration, spacing, length and scaling all agree"))

    print("\nclaims:")
    for text, verdict, why in claims:
        tag = {OK: "[PASS]   ", BAD: "[BLOCKED]", UNK: "[unknown]"}
        print(f"  {tag[verdict]} {text:<38} {why}")

    if a.json:
        print("\n" + json.dumps(dict(systems=systems,
                                     metrics={k: v for k, v in metrics.items()},
                                     issues=[dict(level=l, kind=k, message=msg)
                                             for l, k, msg in issues]),
                                indent=2, ensure_ascii=False, default=str))

    n_blocked = sum(1 for _, v, _ in claims if v == BAD)
    print(f"\n{n_blocked} of {len(claims)} claims blocked")
    return 2 if n_blocked else 0


if __name__ == "__main__":
    sys.exit(main())
