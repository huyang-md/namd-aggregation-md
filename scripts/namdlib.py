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

"""Shared readers for NAMD/CHARMM-GUI systems: PSF, PDB, DCD, XSC.

No external MD packages are required; only numpy.
"""
from __future__ import annotations

import os
import re
import struct

import numpy as np

DT_PS = {"timestep": 2.0}  # fs; overridden per case


# --------------------------------------------------------------------------- #
# PSF
# --------------------------------------------------------------------------- #
def read_psf(path):
    """Return (natoms, [atom dict]) from a CHARMM/XPLOR PSF."""
    with open(path, "r", errors="replace") as fh:
        lines = fh.readlines()
    header = None
    for i, line in enumerate(lines):
        if "!NATOM" in line:
            header = i
            break
    if header is None:
        raise ValueError(f"{path}: no !NATOM block")
    natoms = int(lines[header].split()[0])
    atoms = []
    for line in lines[header + 1: header + 1 + natoms]:
        p = line.split()
        if len(p) < 8:
            continue
        atoms.append(
            dict(
                serial=int(p[0]),
                segid=p[1],
                resid=int(p[2]),
                resname=p[3],
                name=p[4],
                atype=p[5],
                charge=float(p[6]),
                mass=float(p[7]),
            )
        )
    if len(atoms) != natoms:
        raise ValueError(f"{path}: header says {natoms} atoms, parsed {len(atoms)}")
    return natoms, atoms


# --------------------------------------------------------------------------- #
# PDB
# --------------------------------------------------------------------------- #
def read_pdb_box(path):
    """Return (a, b, c, alpha, beta, gamma) from the CRYST1 record, or None."""
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if line.startswith("CRYST1"):
                return (
                    float(line[6:15]), float(line[15:24]), float(line[24:33]),
                    float(line[33:40]), float(line[40:47]), float(line[47:54]),
                )
            if line.startswith(("ATOM", "HETATM")) and not line.startswith("CRYST1"):
                continue
    return None


def read_pdb_atoms(path):
    """Return [(name, resname, resid, x, y, z)] for the first model."""
    out = []
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                out.append((
                    line[12:16].strip(), line[17:20].strip(), int(line[22:26]),
                    float(line[30:38]), float(line[38:46]), float(line[46:54]),
                ))
            elif line.startswith("ENDMDL"):
                break
    return out


# --------------------------------------------------------------------------- #
# NAMD extended system restart (.xsc)
# --------------------------------------------------------------------------- #
def read_xsc(path):
    """Return dict(step=int, cell=3x3 array, origin=array) from a NAMD .xsc file."""
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            p = s.split()
            if len(p) < 7:
                continue
            vals = [float(x) for x in p]
            step = int(vals[0])
            v = np.array(vals[1:])
            cell = np.array([v[0:3], v[3:6], v[6:9]])
            origin = v[9:12]
            return dict(step=step, cell=cell, origin=origin)
    raise ValueError(f"{path}: no data line")


# --------------------------------------------------------------------------- #
# DCD (NAMD / CHARMM binary trajectory)
# --------------------------------------------------------------------------- #
class Dcd:
    """Random-access reader for NAMD DCD trajectories.

    Frame layout (little endian):
        int 48 | 6 x float64 unit cell | int 48
        then, per dimension: int (4*natoms) | natoms x float32 | int (4*natoms)
    """

    def __init__(self, path):
        self.path = path
        self.size = os.path.getsize(path)
        with open(path, "rb") as fh:
            head = fh.read(4)
            if len(head) < 4:
                raise ValueError(f"{path}: empty file")
            first_len = struct.unpack("<i", head)[0]
            if first_len != 84:
                raise ValueError(f"{path}: not a NAMD DCD (record 0 length {first_len})")
            fh.seek(4)
            block = fh.read(84)
            if block[:4] != b"CORD":
                raise ValueError(f"{path}: first record is not CORD ({block[:4]!r})")
            self.nset, self.istart, self.nsavc, self.nstep = struct.unpack("<4i", block[4:20])
            fh.seek(92)
            struct.unpack("<i", fh.read(4))
            ntitle = struct.unpack("<i", fh.read(4))[0]
            fh.seek(100 + ntitle * 80 + 4)
            alen = struct.unpack("<i", fh.read(4))[0]
            self.natoms = struct.unpack("<i", fh.read(4))[0]
            fh.seek(100 + ntitle * 80 + 4 + 4 + 4 + alen)
            self.data_start = fh.tell()
            # Trajectories written by VMD may omit the 48-byte unit-cell record,
            # in which case a frame starts straight at the x-coordinate record.
            probe = fh.read(4)
            self.has_cell = len(probe) == 4 and struct.unpack("<i", probe)[0] == 48
        self.ntitle = ntitle
        self.frame_bytes = ((56 if self.has_cell else 0)
                            + 3 * (8 + 4 * self.natoms))
        avail = self.size - self.data_start
        self.frames = max(0, avail // self.frame_bytes)
        self.partial_bytes = max(0, avail - self.frames * self.frame_bytes)
        self._fh = open(path, "rb")

    # -- helpers ----------------------------------------------------------- #
    def close(self):
        try:
            self._fh.close()
        except Exception:
            pass

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def _read_cell(self, pos):
        """Parse the unit-cell record.

        NAMD writes it as 6 float64 values holding the upper triangle of the cell
        matrix (h11 h12 h22 h13 h23 h33); CHARMM-style writers may use 9. Both are
        handled. Orthorhombic boxes come out diagonal, which is the common case.
        """
        if not self.has_cell:
            return None
        fh = self._fh
        fh.seek(pos)
        reclen = struct.unpack("<i", fh.read(4))[0]
        raw = fh.read(reclen)
        fh.read(4)
        if reclen == 72:
            cell = np.frombuffer(raw, dtype="<f8").reshape(3, 3).copy()
        elif reclen == 48:
            v = np.frombuffer(raw, dtype="<f8")
            cell = np.array([[v[0], v[1], v[3]],
                             [0.0, v[2], v[4]],
                             [0.0, 0.0, v[5]]], dtype=float)
        else:
            raise ValueError(f"{self.path}: unexpected unit-cell record length {reclen}")
        return cell

    def cell(self, i):
        return self._read_cell(self.data_start + i * self.frame_bytes)

    def frame(self, i, indices=None):
        """Return (coords[n,3] float32, cell 3x3 or None) for frame i."""
        if i < 0:
            i += self.frames
        if not (0 <= i < self.frames):
            raise IndexError(f"frame {i} out of range (0..{self.frames - 1})")
        pos = self.data_start + i * self.frame_bytes
        cell = self._read_cell(pos)
        n = self.natoms
        out = np.empty((len(indices) if indices is not None else n, 3), dtype=np.float32)
        fh = self._fh
        base = pos + (56 if self.has_cell else 0)
        for d in range(3):
            fh.seek(base + d * (8 + 4 * n) + 4)
            block = np.frombuffer(fh.read(4 * n), dtype="<f4")
            out[:, d] = block if indices is None else block[indices]
        return out, cell

    def frames_range(self, indices=None, stride=1, start=0, stop=None):
        stop = self.frames if stop is None else min(stop, self.frames)
        for i in range(start, stop, max(1, stride)):
            yield i, *self.frame(i, indices)


# --------------------------------------------------------------------------- #
# Geometry / bookkeeping
# --------------------------------------------------------------------------- #
def frame_bytes_for(natoms):
    return 56 + 3 * (8 + 4 * natoms)


def trajectory_bytes(natoms, ns, ps_per_frame):
    frames = int(round(ns * 1000.0 / ps_per_frame))
    return frames * frame_bytes_for(natoms)


def human_bytes(n):
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if abs(n) < 1024.0 or unit == "TB":
            return f"{n:,.1f} {unit}"
        n /= 1024.0


def minimum_image(delta, cell):
    """Wrap displacement(s) into [-L/2, L/2) for an orthorhombic cell."""
    lengths = np.linalg.norm(cell, axis=1)
    return delta - np.round(delta / lengths) * lengths


def unwrap_molecule(x, cell):
    """Return the molecule's atoms in the global frame, made whole across the boundary.

    Atoms are moved by whole box vectors so that they sit within half a box of the
    molecule's first atom. The result stays where the molecule actually is: it must
    not be re-centred on the first atom, or molecules can no longer be compared to
    each other.
    """
    lengths = np.linalg.norm(cell, axis=1)
    ref = x[0]
    d = x - ref
    d -= np.round(d / lengths) * lengths
    return d + ref


def plane(pts):
    """Best-fit plane: return (centroid, unit normal, rms deviation)."""
    c = pts.mean(axis=0)
    _, _, vt = np.linalg.svd(pts - c, full_matrices=False)
    n = vt[2]
    rms = float(np.sqrt((((pts - c) @ n) ** 2).mean()))
    return c, n, rms


NAME_RE = re.compile(r"^[A-Za-z0-9_+\-']+$")


def parse_name_list(text):
    return [t.strip() for t in re.split(r"[,\s]+", text.strip()) if t.strip()]
