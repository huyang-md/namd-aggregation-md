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

"""Zero-dependency demo: synthesise a small self-assembly case, then analyse it.

Builds examples/demo-case/ (PSF + PDB + DCD, solute only): 40 dye molecules
(8 atoms each) in a 100 A box -- 20 of them pre-assembled into a spherical
cluster at the centre, 20 dispersed freely -- then runs scripts/aggregate.py
on it exactly the way a user would, so what you see is the real tool, not a
canned transcript.

    python examples/demo.py            # from the repository root

Requires only Python 3.10+ and numpy -- no NAMD, VMD or Packmol.

Two ways to play with it:

* edit the constants below (cluster size, box, jitter, frames) and re-run;
* re-analyse the generated case with different criteria:

    python scripts/aggregate.py --dir examples/demo-case --traj 04_prod.dcd --mode centroid
    python scripts/aggregate.py --dir examples/demo-case --traj 04_prod.dcd --cutoff 6

Outputs land in examples/demo-case/analysis/; start with summary.txt. The run
also writes clusters_vmd.tcl (--vmd-out): load ionized.psf + 04_prod.dcd in
VMD and source that file to see the molecules coloured by cluster.
"""
from __future__ import annotations

import os
import struct
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CASE = os.path.join(HERE, "demo-case")

BOX = 100.0            # cubic box edge in A
CENTRE = np.array([BOX / 2.0] * 3)   # the pre-assembled ball sits at the box centre
N_CLUSTER = 20         # molecules pre-assembled into a ball at the box centre
N_FREE = 20            # molecules dispersed in the box
CLUSTER_RADIUS = 8.0   # radius of the pre-assembled ball in A
MIN_SEP = 4.8          # minimum centre-centre distance inside the ball
ATOMS_PER_MOL = 8
CUBE_HALF = 1.25       # a molecule is 8 atoms on a cube of side 2.5 A
N_FRAMES = 21
NSAVC = 100000         # steps between frames -> 200 ps per frame (2 fs timestep)
JITTER = 0.25          # thermal jitter, sigma in A
SEED = 20260101


# --------------------------------------------------------------------------- #
# geometry
# --------------------------------------------------------------------------- #
def _rejection_sample(rng, n, try_make):
    """Call try_make(pts) until it has produced n points (None = rejected draw)."""
    pts = []
    for _ in range(200000):
        p = try_make(pts)
        if p is not None:
            pts.append(p)
            if len(pts) == n:
                return np.array(pts)
    raise RuntimeError("rejection sampling stalled; lower MIN_SEP and re-run")


def _cluster_centers(rng):
    sep = MIN_SEP

    def try_make(pts):
        nonlocal sep
        p = CENTRE + rng.uniform(-CLUSTER_RADIUS, CLUSTER_RADIUS, 3)
        if np.linalg.norm(p - CENTRE) > CLUSTER_RADIUS:
            return None
        if pts and min(np.linalg.norm(p - q) for q in pts) < sep:
            return None
        return p

    try:
        return _rejection_sample(rng, N_CLUSTER, try_make)
    except RuntimeError:
        sep = 4.0     # a denser ball still clusters fine at the 4 A contact cutoff
        return _rejection_sample(rng, N_CLUSTER, try_make)


def _free_centers(rng):
    def try_make(pts):
        p = rng.uniform(10.0, BOX - 10.0, 3)
        if np.linalg.norm(p - CENTRE) < 30.0:
            return None
        if pts and min(np.linalg.norm(p - q) for q in pts) < 14.0:
            return None
        return p

    return _rejection_sample(rng, N_FREE, try_make)


def _molecule_atoms(center):
    offs = np.array([(sx * CUBE_HALF, sy * CUBE_HALF, sz * CUBE_HALF)
                     for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)])
    return center + offs


# --------------------------------------------------------------------------- #
# synthetic file writers (same on-disk formats as tests/test_core.py)
# --------------------------------------------------------------------------- #
def _write_psf(path, nmol):
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write("PSF\n\n       1 !NTITLE\n REMARKS synthetic demo case\n\n")
        fh.write(f"{nmol * ATOMS_PER_MOL:8d} !NATOM\n")
        serial = 0
        for resid in range(1, nmol + 1):
            for k in range(ATOMS_PER_MOL):
                serial += 1
                name = f"C{k + 1}"
                fh.write(f"{serial:8d} S       {resid:4d} DYE  {name:<4s} C    "
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
            resid = (i - 1) // ATOMS_PER_MOL + 1
            name = f"C{(i - 1) % ATOMS_PER_MOL + 1}"
            fh.write("ATOM  %5d %-4s %-3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00\n"
                     % (i, name, "DYE", resid, x, y, z))
        fh.write("END\n")


def _write_dcd(path, frames):
    nframes, natoms = len(frames), len(frames[0])
    title = b"synthetic demo".ljust(80)
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


# --------------------------------------------------------------------------- #
def main():
    rng = np.random.default_rng(SEED)
    os.makedirs(CASE, exist_ok=True)

    centers = np.vstack([_cluster_centers(rng), _free_centers(rng)])
    nmol = len(centers)
    base = np.vstack([_molecule_atoms(c) for c in centers])
    print(f"case      : {os.path.relpath(CASE, ROOT)}", flush=True)
    print(f"molecules : {nmol} x DYE ({ATOMS_PER_MOL} atoms) = {nmol * ATOMS_PER_MOL} atoms "
          f"in a {BOX:.0f} A box", flush=True)
    print(f"setup     : {N_CLUSTER} pre-assembled (ball r={CLUSTER_RADIUS:.0f} A at the centre) "
          f"+ {N_FREE} free", flush=True)
    print(f"trajectory: {N_FRAMES} frames {NSAVC * 2 / 1000:.0f} ps apart "
          f"= {N_FRAMES * NSAVC * 2 / 1000 / 1000:.1f} ns\n", flush=True)

    _write_psf(os.path.join(CASE, "ionized.psf"), nmol)
    _write_pdb(os.path.join(CASE, "ionized.pdb"), base)
    frames = [base + rng.normal(0.0, JITTER, base.shape) for _ in range(N_FRAMES)]
    _write_dcd(os.path.join(CASE, "04_prod.dcd"), frames)

    print("running scripts/aggregate.py on the case:\n", flush=True)
    script = os.path.join(ROOT, "scripts", "aggregate.py")
    r = subprocess.run([sys.executable, script, "--dir", CASE, "--traj", "04_prod.dcd",
                        "--vmd-out"])
    if r.returncode != 0:
        sys.exit(r.returncode)

    print("\nnext steps -- play with the analysis on this case:")
    agg = os.path.join(ROOT, "scripts", "aggregate.py")
    print(f'  python "{agg}" --dir "{CASE}" --traj 04_prod.dcd --mode centroid --cutoff 6')
    print(f'  python "{agg}" --dir "{CASE}" --traj 04_prod.dcd --cutoff 6')
    print("  see the aggregates in VMD: load ionized.psf + 04_prod.dcd, then source")
    print(f'  "{os.path.join(CASE, "analysis", "clusters_vmd.tcl")}"')
    print(f'  outputs are in {os.path.join(CASE, "analysis")} -- start with summary.txt')
    return 0


if __name__ == "__main__":
    sys.exit(main())
