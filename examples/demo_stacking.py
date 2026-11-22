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

"""Zero-dependency demo for the stacking-geometry layer (analyze.py).

Builds examples/demo-stacking-case/: 12 'LIG' molecules, each a flat six-atom
hexagon, arranged as two columns of six -- 3.5 A between neighbouring planes,
4 A lateral slip between neighbours. That is the textbook geometry of a
J-aggregate, and analyze.py should find exactly that: slip angle around
49 deg (J-type), interplanar distance around 3.5 A, nearest molecular
centres around 5.3 A. Then it runs scripts/analyze.py the way a user would,
so what you see is the real tool, not a canned transcript.

    python examples/demo_stacking.py   # from the repository root

Requires only Python 3.10+ and numpy -- no NAMD, VMD or Packmol.

analyze.py reads its atom groups from case.json, so the demo also writes a
minimal one (plane_groups / plane_pair / rdf_pairs) -- the same fields that
inspect_molecule.py fills in on a real case. Play with it:

    python scripts/analyze.py --dir examples/demo-stacking-case --traj 04_prod.dcd --slip-threshold 60
    python scripts/analyze.py --dir examples/demo-stacking-case --traj 04_prod.dcd --no-rdf

Outputs land in examples/demo-stacking-case/analysis/; start with summary.txt.
"""
from __future__ import annotations

import json
import os
import struct
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CASE = os.path.join(HERE, "demo-stacking-case")

BOX = 100.0            # cubic box edge in A
RING_ATOMS = 6         # one flat six-atom hexagon per molecule
RING_RADIUS = 1.4      # hexagon radius in A
N_STACKS = 2
PER_STACK = 6
D_PLANE = 3.5          # interplanar spacing along z (A)
SLIP_XY = 4.0          # lateral slip between neighbours along x (A)
STACK1_XYZ = (30.0, 50.0, 40.0)
STACK2_XYZ = (70.0, 50.0, 40.0)
N_FRAMES = 21
NSAVC = 100000         # steps between frames -> 200 ps per frame (2 fs timestep)
JITTER = 0.08          # thermal jitter, sigma in A
SEED = 20260101


# --------------------------------------------------------------------------- #
# geometry
# --------------------------------------------------------------------------- #
def _molecule_atoms(center, slip_steps):
    """One flat hexagon: atoms at RING_RADIUS in xy, z = 0 inside the molecule.

    slip_steps shifts the molecule along x, which is what makes neighbours
    J-type stacked instead of face-to-face.
    """
    angles = np.arange(RING_ATOMS) * (2.0 * np.pi / RING_ATOMS)
    local = np.column_stack([RING_RADIUS * np.cos(angles),
                             RING_RADIUS * np.sin(angles),
                             np.zeros(RING_ATOMS)])
    return np.asarray(center) + np.array([slip_steps * SLIP_XY, 0.0, 0.0]) + local


def _all_centers():
    centers = []
    for base in (STACK1_XYZ, STACK2_XYZ):
        for n in range(PER_STACK):
            centers.append((np.asarray(base) + np.array([0.0, 0.0, n * D_PLANE]), n))
    return centers


# --------------------------------------------------------------------------- #
# synthetic file writers (same on-disk formats as tests/test_core.py)
# --------------------------------------------------------------------------- #
def _write_psf(path, nmol):
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write("PSF\n\n       1 !NTITLE\n REMARKS synthetic stacking demo\n\n")
        fh.write(f"{nmol * RING_ATOMS:8d} !NATOM\n")
        serial = 0
        for resid in range(1, nmol + 1):
            for k in range(RING_ATOMS):
                serial += 1
                name = f"C{k + 1}"
                fh.write(f"{serial:8d} S       {resid:4d} LIG  {name:<4s} C    "
                         f"{0.0:10.6f} {12.010:10.4f}           0\n")
        fh.write("\n       0 !NBOND\n")
        fh.write("\n       0 !NIMPHI: impropers\n")
        fh.write("\n       0 !NDON: donors\n")
        fh.write("\n       0 !NACC: acceptors\n")
        fh.write("\n       0 !NNB\n")


def _write_pdb(path, coords):
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n"
                 % (BOX, BOX, BOX, 90.0, 90.0, 90.0))
        for i, (x, y, z) in enumerate(coords, 1):
            resid = (i - 1) // RING_ATOMS + 1
            name = f"C{(i - 1) % RING_ATOMS + 1}"
            fh.write("ATOM  %5d %-4s %-3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00\n"
                     % (i, name, "LIG", resid, x, y, z))
        fh.write("END\n")


def _write_dcd(path, frames):
    nframes, natoms = len(frames), len(frames[0])
    title = b"synthetic stacking demo".ljust(80)
    buf = bytearray()
    buf += struct.pack("<i", 84)
    buf += b"CORD" + struct.pack("<4i", nframes, 0, NSAVC, nframes * NSAVC) + b"\0" * (84 - 4 - 16)
    buf += struct.pack("<i", 84)
    buf += struct.pack("<i", 84) + struct.pack("<i", 1) + title + struct.pack("<i", 84)
    buf += struct.pack("<i", 4) + struct.pack("<i", natoms) + struct.pack("<i", 4)
    for coords in frames:
        buf += struct.pack("<i", 48)
        # NAMD order: h11 h12 h22 h13 h23 h33
        buf += struct.pack("<6d", BOX, 0.0, BOX, 0.0, 0.0, BOX)
        buf += struct.pack("<i", 48)
        for d in range(3):
            data = np.asarray(coords, dtype="<f4")[:, d].tobytes()
            buf += struct.pack("<i", len(data)) + data + struct.pack("<i", len(data))
    with open(path, "wb") as fh:
        fh.write(bytes(buf))


def _write_case(path):
    case = {
        "name": "demo-stacking",
        "psf": "ionized.psf",
        "pdb": "ionized.pdb",
        "ligand_resname": "LIG",
        "timestep_fs": 2.0,
        "analysis": {
            "plane_groups": {"ring1": [f"C{k + 1}" for k in range(RING_ATOMS)]},
            "plane_pair": ["ring1", "ring1"],
            "rdf_pairs": [["ring:ring1", "ring:ring1"]],
            "cluster_cutoff_A": 4.0,
        },
    }
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(case, fh, indent=2)
        fh.write("\n")


# --------------------------------------------------------------------------- #
def main():
    rng = np.random.default_rng(SEED)
    os.makedirs(CASE, exist_ok=True)

    centers = _all_centers()
    nmol = len(centers)
    base = np.vstack([_molecule_atoms(c, slip) for c, slip in centers])
    print(f"case      : {os.path.relpath(CASE, ROOT)}")
    print(f"molecules : {nmol} x LIG ({RING_ATOMS}-atom hexagons) = {nmol * RING_ATOMS} atoms "
          f"in a {BOX:.0f} A box", flush=True)
    print(f"setup     : {N_STACKS} columns of {PER_STACK}, planes {D_PLANE} A apart, "
          f"neighbours slipped {SLIP_XY} A sideways (J-aggregate geometry)", flush=True)
    print(f"trajectory: {N_FRAMES} frames {NSAVC * 2 / 1000:.0f} ps apart "
          f"= {N_FRAMES * NSAVC * 2 / 1000 / 1000:.1f} ns", flush=True)
    print("expected  : slip angle about 49 deg (J-type), interplanar d about 3.5 A, "
          "nearest centres about 5.3 A\n", flush=True)

    _write_psf(os.path.join(CASE, "ionized.psf"), nmol)
    _write_pdb(os.path.join(CASE, "ionized.pdb"), base)
    _write_case(os.path.join(CASE, "case.json"))
    frames = [base + rng.normal(0.0, JITTER, base.shape) for _ in range(N_FRAMES)]
    _write_dcd(os.path.join(CASE, "04_prod.dcd"), frames)

    print("running scripts/analyze.py on the case:\n", flush=True)
    script = os.path.join(ROOT, "scripts", "analyze.py")
    r = subprocess.run([sys.executable, script, "--dir", CASE, "--traj", "04_prod.dcd"])
    if r.returncode != 0:
        sys.exit(r.returncode)

    print("\nnext steps -- play with the analysis on this case:")
    agg = os.path.join(ROOT, "scripts", "analyze.py")
    print(f'  python "{agg}" --dir "{CASE}" --traj 04_prod.dcd --slip-threshold 60')
    print("      (60 deg turns the same geometry 'H-type' -- thresholds decide labels,")
    print("       which is why you must state your threshold in papers)")
    print(f'  python "{agg}" --dir "{CASE}" --traj 04_prod.dcd --no-rdf')
    print("  outputs are in " + os.path.join(CASE, "analysis") + " -- start with summary.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
