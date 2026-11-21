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

"""Compare sibling systems for methodological comparability.

For a paper that compares aggregation across substituents, the systems must agree on
things that are easy to get wrong when each one was assembled by hand: solute
concentration (Packmol box + VMD solvate padding decide it, not the intent), frame
spacing, production length, and the force-field knobs that silently change physics.

    python compare_systems.py --root <parent dir>
    python compare_systems.py --root <parent dir> --resname LIG --json
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
import preflight  # noqa: E402
import scaffold  # noqa: E402

PROD_RE = re.compile(r"prod", re.I)
DERIVED_RE = re.compile(r"nowater|no_?water|dry|solute|ligand_?only|small|stripped|\.bak$", re.I)
AVOGADRO_PER_A3_PER_M = 6.02214076e-4  # (mol/L) -> molecules per A^3


def production_trajectories(directory):
    """Real MD output only: water-stripped VMD exports are not production runs."""
    out = []
    for f in glob.glob(os.path.join(directory, "*.dcd")):
        base = os.path.basename(f)
        if DERIVED_RE.search(base) or not PROD_RE.search(base):
            continue
        out.append(f)
    return out


def analyse_case(directory, resname):
    info = dict(name=os.path.basename(directory.rstrip("\\/")), dir=directory)

    case = None
    cpath = os.path.join(directory, scaffold.CASE_FILE)
    if os.path.exists(cpath):
        try:
            case = scaffold.load_case(directory)
        except Exception:  # noqa: BLE001
            case = None
    psf_name = (case or {}).get("psf", "ionized.psf")
    pdb_name = (case or {}).get("pdb", "ionized.pdb")
    psf = os.path.join(directory, psf_name)
    pdb = os.path.join(directory, pdb_name)
    if not os.path.exists(psf):
        return None

    natoms, atoms = namdlib.read_psf(psf)
    info["natoms"] = natoms
    lig = [a for a in atoms if a["resname"] == resname]
    if lig:
        nmol = len({a["resid"] for a in lig})
        info["molecules"] = nmol
        info["atoms_per_molecule"] = len(lig) // nmol if nmol else 0

    # production span and the box actually used during production
    best = None
    spacings = set()
    vol_prod = None
    for f in production_trajectories(directory):
        try:
            d = namdlib.Dcd(f)
        except Exception:  # noqa: BLE001
            continue
        nsavc = d.nsavc
        if nsavc:
            spacings.add(nsavc * 2.0 / 1000.0)
        end_time = ((d.istart + d.frames * nsavc) * 2.0 / 1e6) if d.frames else 0.0
        if best is None or end_time > best[0]:
            cell = None
            try:
                cell = d.cell(d.frames - 1) if d.frames else None
            except Exception:  # noqa: BLE001
                cell = None
            best = [end_time, cell]
        d.close()
    if best:
        info["prod_ns"] = round(best[0], 1)
        if best[1] is not None:
            vol_prod = abs(float(np_linalg_det(best[1])))
    if spacings:
        info["spacing_ps"] = sorted(spacings)

    cryst = namdlib.read_pdb_box(pdb) if os.path.exists(pdb) else None
    vol_init = cryst[0] * cryst[1] * cryst[2] if cryst else None
    if cryst:
        info["box"] = [round(x, 2) for x in cryst[:3]]
        info["volume_initial_nm3"] = round(vol_init / 1000.0, 1)
    vol = vol_prod or vol_init
    if vol:
        info["volume_nm3"] = round(vol / 1000.0, 1)
        info["volume_source"] = "last production frame" if vol_prod else "initial CRYST1"
        if info.get("molecules"):
            info["conc_mM"] = round(info["molecules"] / (vol * AVOGADRO_PER_A3_PER_M) * 1000.0, 1)

    # physics knobs, per .conf; production and equilibration may disagree
    variants = {}
    prod_conf = None
    best_run = -1.0
    for f in sorted(glob.glob(os.path.join(directory, "*.conf"))):
        try:
            scalars, _ = preflight.parse_conf(f)
        except Exception:  # noqa: BLE001
            continue
        if "oneFourScaling" in scalars:
            try:
                variants[os.path.basename(f)] = float(scalars["oneFourScaling"])
            except ValueError:
                pass
        try:
            run = float(scalars.get("run", "0"))
        except ValueError:
            run = 0.0
        base = os.path.basename(f).lower()
        if ("prod" in base or "product" in base) and run >= best_run:
            best_run = run
            prod_conf = scalars
    if prod_conf:
        for key, label in (("oneFourScaling", "scaling_14"), ("timestep", "timestep_fs"),
                           ("cutoff", "cutoff")):
            if key in prod_conf:
                try:
                    info[label] = float(prod_conf[key])
                except ValueError:
                    pass
    if variants:
        info["scaling_by_conf"] = variants
    return info


def np_linalg_det(m):
    (a, b, c), (d, e, f), (g, h, i) = m
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def compare(rows):
    issues = []
    concs = [(r["name"], r["conc_mM"]) for r in rows if r.get("conc_mM")]
    if len(concs) > 1:
        lo = min(concs, key=lambda t: t[1])
        hi = max(concs, key=lambda t: t[1])
        ratio = hi[1] / lo[1] if lo[1] else 0
        if ratio > 1.05:
            issues.append(("FAIL", "concentration",
                           f"{hi[0]} is {ratio:.2f}x more concentrated than {lo[0]} "
                           f"({hi[1]} vs {lo[1]} mM) - aggregate size differences may follow "
                           "the concentration, not the substituent"))
        else:
            issues.append(("INFO", "concentration",
                           f"concentrations within {ratio:.2f}x"))

    sp = {r["name"]: tuple(r["spacing_ps"]) for r in rows if r.get("spacing_ps")}
    if len(set(sp.values())) > 1:
        detail = ", ".join(f"{k}={'/'.join(str(x) for x in v)}" for k, v in sorted(sp.items()))
        issues.append(("FAIL", "spacing", f"frame spacing differs across systems: {detail}"))

    pl = {r["name"]: r["prod_ns"] for r in rows if r.get("prod_ns")}
    if pl:
        lo = min(pl.values())
        hi = max(pl.values())
        if hi > 0 and (hi - lo) / hi > 0.02:
            issues.append(("FAIL", "length",
                           f"production length differs: {lo} to {hi} ns"))
        elif hi:
            issues.append(("INFO", "length", f"production length {lo}-{hi} ns"))

    for label, pretty in (("scaling_14", "1-4scaling"), ("timestep_fs", "timestep"),
                          ("cutoff", "cutoff")):
        vals = {r["name"]: r[label] for r in rows if label in r}
        if len(set(vals.values())) > 1:
            detail = ", ".join(f"{k}={v:g}" for k, v in sorted(vals.items()))
            issues.append(("FAIL", label, f"{pretty} differs: {detail}"))

    for r in rows:
        by_conf = r.get("scaling_by_conf", {})
        if len({round(v, 6) for v in by_conf.values()}) > 1:
            detail = ", ".join(f"{k}={v:g}" for k, v in sorted(by_conf.items()))
            issues.append(("FAIL", "scaling-within",
                           f"{r['name']}: 1-4scaling changes between its own stages -> {detail}"))

    nmol = {r["name"]: r.get("molecules") for r in rows if r.get("molecules")}
    if len(set(nmol.values())) > 1:
        detail = ", ".join(f"{k}={v}" for k, v in sorted(nmol.items()))
        issues.append(("FAIL", "molecules", f"molecule count differs: {detail}"))
    return issues


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--root", required=True, help="parent directory holding one system per subdirectory")
    p.add_argument("--resname", default="LIG")
    p.add_argument("--json", action="store_true")
    a = p.parse_args(argv)

    root = os.path.abspath(a.root)
    rows = []
    for name in sorted(os.listdir(root)):
        d = os.path.join(root, name)
        if not os.path.isdir(d):
            continue
        try:
            r = analyse_case(d, a.resname)
        except Exception as exc:  # noqa: BLE001
            print(f"[skip] {name}: {exc}")
            continue
        if r:
            rows.append(r)

    if not rows:
        sys.exit(f"no systems with a PSF found under {root}")
    issues = compare(rows)

    if a.json:
        print(json.dumps(dict(root=root, systems=rows,
                              issues=[dict(level=l, kind=k, message=m) for l, k, m in issues]),
                         indent=2, ensure_ascii=False))
        return 2 if any(l == "FAIL" for l, _, _ in issues) else 0

    hdr = (f"{'system':<14}{'molecules':>10}{'atoms':>10}{'box (A)':>26}"
           f"{'V (nm3)':>10}{'conc (mM)':>11}{'prod (ns)':>11}{'ps/frame':>10}")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        box = " x ".join(f"{x:.1f}" for x in r.get("box", [])) or "-"
        sp = "/".join(f"{x:g}" for x in r.get("spacing_ps", [])) or "-"
        print(f"{r['name']:<14}{r.get('molecules', '-'):>10}{r.get('natoms', '-'):>10,}"
              f"{box:>26}{r.get('volume_nm3', '-'):>10}{r.get('conc_mM', '-'):>11}"
              f"{r.get('prod_ns', '-'):>11}{sp:>10}")

    print("\ncomparability:")
    for level, kind, message in issues:
        tag = {"FAIL": "[FAIL]", "WARN": "[warn]", "INFO": "[ok]  "}[level]
        print(f"  {tag} {message}")
    n_fail = sum(1 for l, _, _ in issues if l == "FAIL")
    print(f"\n{n_fail} blocking difference(s) across {len(rows)} systems")
    return 2 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
