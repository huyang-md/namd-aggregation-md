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

"""Discover the aromatic ring systems of a solute and propose analysis settings.

Any small-molecule aggregation case can be analysed without hard-coding atom names:
this reads the topology, finds rings, groups fused rings into ring systems, checks
planarity from the coordinates, and writes the plane groups / RDF sites into case.json.

    python inspect_molecule.py --dir <case dir>
    python inspect_molecule.py --dir <case dir> --write
    python inspect_molecule.py --dir <case dir> --resname LIG --min-rms 0.2
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import namdlib  # noqa: E402
import scaffold  # noqa: E402

WATER = {"TIP3", "TIP2", "WAT", "HOH", "SOL", "TIPS", "SPC"}
IONS = {"SOD", "CLA", "POT", "CAL", "MG", "ZN", "NA", "CL", "K", "CA", "CES", "BRO", "IOD"}
MAX_RING = 7


# --------------------------------------------------------------------------- #
def read_psf_bonds(path):
    lines = open(path, "r", errors="replace").read().splitlines()
    nb_i = next((i for i, l in enumerate(lines) if "!NBOND" in l), None)
    if nb_i is None:
        return []
    nbonds = int(lines[nb_i].split()[0])
    flat = []
    for ln in lines[nb_i + 1:]:
        flat.extend(ln.split())
        if len(flat) >= 2 * nbonds:
            break
    flat = flat[:2 * nbonds]
    return [(int(flat[k]) - 1, int(flat[k + 1]) - 1) for k in range(0, len(flat), 2)]


def find_rings(adj, natoms, maxlen=MAX_RING):
    """All simple cycles up to maxlen, deduplicated by atom set."""
    rings = set()

    def walk(start, cur, path, seen):
        if len(path) > maxlen:
            return
        for nb in adj[cur]:
            if nb == start and len(path) >= 3:
                rings.add(frozenset(path))
            elif nb not in seen and nb > start:
                seen.add(nb)
                walk(start, nb, path + [nb], seen)
                seen.discard(nb)

    for s in range(natoms):
        walk(s, s, [s], {s})
    return [sorted(r) for r in rings]


def ring_systems(rings):
    """Fuse rings that share an edge (>=2 atoms) into systems."""
    parent = list(range(len(rings)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(len(rings)):
        for j in range(i + 1, len(rings)):
            if len(set(rings[i]) & set(rings[j])) >= 2:
                ri, rj = find(i), find(j)
                if ri != rj:
                    parent[rj] = ri

    groups = defaultdict(list)
    for i, r in enumerate(rings):
        groups[find(i)].extend(r)
    return [sorted(set(v)) for v in groups.values()]


def detect_solute(atoms, resname=None):
    if resname:
        picked = [a for a in atoms if a["resname"] == resname]
        return resname, picked
    counts = defaultdict(int)
    for a in atoms:
        rn = a["resname"]
        if rn in WATER or rn in IONS:
            continue
        counts[rn] += 1
    if not counts:
        sys.exit("no solute residue found (everything looks like water/ions); pass --resname")
    best = max(counts, key=lambda k: counts[k])
    return best, [a for a in atoms if a["resname"] == best]


def planarity(coords):
    c = coords.mean(axis=0)
    _, _, vt = np.linalg.svd(coords - c, full_matrices=False)
    n = vt[2]
    return float(np.sqrt((((coords - c) @ n) ** 2).mean()))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dir", required=True)
    p.add_argument("--psf")
    p.add_argument("--pdb")
    p.add_argument("--resname", help="solute residue name (default: auto-detect)")
    p.add_argument("--min-ring", type=int, default=5, help="smallest ring to consider")
    p.add_argument("--max-rms", type=float, default=0.2,
                   help="planarity threshold (A) for a ring system to be usable as a plane")
    p.add_argument("--write", action="store_true", help="merge the proposal into case.json")
    p.add_argument("--json", action="store_true")
    a = p.parse_args(argv)

    d = os.path.abspath(a.dir)
    case = scaffold.load_case(d) if os.path.exists(os.path.join(d, scaffold.CASE_FILE)) else None
    psf = a.psf or os.path.join(d, (case or {}).get("psf", "ionized.psf"))
    pdb = a.pdb or os.path.join(d, (case or {}).get("pdb", "ionized.pdb"))
    if not os.path.exists(psf):
        sys.exit(f"missing {psf}")

    natoms, atoms = namdlib.read_psf(psf)
    resname, solute = detect_solute(atoms, a.resname or (case or {}).get("ligand_resname"))
    copies = len({x["resid"] for x in solute})
    per = len(solute) // copies if copies else 0
    print(f"psf       : {psf}")
    print(f"solute    : {resname}  ({copies} copies, {per} atoms each) of {natoms:,} total")
    if per == 0:
        sys.exit("solute residue has no atoms")

    first = [x for x in solute if x["resid"] == sorted({x["resid"] for x in solute})[0]]
    local = {x["serial"] - 1: k for k, x in enumerate(first)}
    bonds = [(local[i], local[j]) for i, j in read_psf_bonds(psf) if i in local and j in local]
    adj = defaultdict(set)
    for i, j in bonds:
        adj[i].add(j)
        adj[j].add(i)

    rings = [r for r in find_rings(adj, per) if len(r) >= a.min_ring]
    if not rings:
        sys.exit("no rings found — is the solute cyclic? (linear molecules need explicit groups)")
    systems = ring_systems(rings)

    pdb_atoms = namdlib.read_pdb_atoms(os.path.join(pdb) if os.path.exists(pdb) else psf)
    if len(pdb_atoms) < per:
        sys.exit(f"{pdb} has {len(pdb_atoms)} atoms, expected at least {per}")
    xyz = np.array([[t[3], t[4], t[5]] for t in pdb_atoms[:per]])

    names = [x["name"] for x in first]
    types = [x["atype"] for x in first]

    proposal = {}
    rows = []
    for k, sys_atoms in enumerate(sorted(systems, key=lambda s: (-len(s), s[0])), 1):
        label = f"ring{k}"
        arom = sum(1 for i in sys_atoms if "R" in types[i])
        rms = planarity(xyz[sys_atoms])
        ok = rms <= a.max_rms
        rows.append((label, len(sys_atoms), f"{arom}/{len(sys_atoms)}", rms, ok,
                     [names[i] for i in sys_atoms]))
        if ok:
            proposal[label] = [names[i] for i in sys_atoms]

    hdr = f"{'group':<8}{'atoms':>7}{'aromatic':>10}{'planarity rms':>15}{'usable':>8}  atoms"
    print("\n" + hdr)
    print("-" * len(hdr))
    for label, n, arom, rms, ok, ats in rows:
        shown = ",".join(ats[:8]) + (" ..." if len(ats) > 8 else "")
        print(f"{label:<8}{n:>7}{arom:>10}{rms:>15.3f}{('yes' if ok else 'no'):>8}  {shown}")

    usable = [label for label, _, _, _, ok, _ in rows if ok]
    if not usable:
        print(f"\nno ring system is planar within {a.max_rms} A; raise --max-rms or "
              "define plane_groups by hand in case.json")
    else:
        primary = max(usable, key=lambda l: len(proposal[l]))
        analysis = {
            "plane_groups": proposal,
            "plane_pair": [primary, primary],
            "rdf_pairs": [[f"ring:{primary}", f"ring:{primary}"]],
        }
        print(f"\nproposal (primary plane: {primary})")
        print(json.dumps(analysis, indent=2, ensure_ascii=False))

        if a.write:
            if case is None:
                sys.exit("no case.json here; run scaffold.py first or omit --write")
            case["analysis"].update(analysis)
            path = os.path.join(d, scaffold.CASE_FILE)
            with open(path, "w", encoding="utf-8", newline="\n") as fh:
                json.dump(case, fh, indent=2, ensure_ascii=False)
                fh.write("\n")
            print(f"\nwrote {path}")

    if a.json:
        print("\n" + json.dumps(dict(
            solute=resname, copies=copies, atoms_per_copy=per,
            groups=[dict(name=l, atoms=ats, planarity_rms=r,
                         aromatic=ar, usable=ok) for l, _, ar, r, ok, ats in rows]),
            indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())