"""Readers: PSF, PDB, XSC, DCD. Builds synthetic files so no MD data is needed."""
import os
import struct
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
import namdlib  # noqa: E402


# --------------------------------------------------------------------------- #
# synthetic file builders
# --------------------------------------------------------------------------- #
def write_dcd(path, coords_list, istart=0, nsavc=1000, title="synthetic"):
    """Write a NAMD-style DCD. coords_list is a list of (natoms,3) arrays."""
    nframes = len(coords_list)
    natoms = len(coords_list[0])
    titles = title.encode("ascii").ljust(80)
    tl = 4 + 80
    buf = bytearray()
    buf += struct.pack("<i", 84)
    buf += b"CORD" + struct.pack("<4i", nframes, istart, nsavc, istart + nframes * nsavc) + b"\0" * (84 - 4 - 16)
    buf += struct.pack("<i", 84)
    buf += struct.pack("<i", tl) + struct.pack("<i", 1) + titles + struct.pack("<i", tl)
    buf += struct.pack("<i", 4) + struct.pack("<i", natoms) + struct.pack("<i", 4)
    for coords in coords_list:
        buf += struct.pack("<i", 48)
        # NAMD order: h11 h12 h22 h13 h23 h33
        buf += struct.pack("<6d", 100.0, 0.0, 100.0, 0.0, 0.0, 100.0)
        buf += struct.pack("<i", 48)
        for d in range(3):
            data = np.asarray(coords, dtype="<f4")[:, d].tobytes()
            buf += struct.pack("<i", len(data)) + data + struct.pack("<i", len(data))
    with open(path, "wb") as fh:
        fh.write(bytes(buf))
    return natoms


def write_psf(path, atoms, bonds=()):
    """atoms: list of (resid, resname, name, atype, charge, mass)."""
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write("PSF\n\n       1 !NTITLE\n REMARKS synthetic\n\n")
        fh.write(f"{len(atoms):8d} !NATOM\n")
        for i, (resid, resname, name, atype, q, m) in enumerate(atoms, 1):
            fh.write(f"{i:8d} S       {resid:4d} {resname:<4s} {name:<4s} {atype:<4s} "
                     f"{q:10.6f} {m:10.4f}           0\n")
        fh.write(f"\n{len(bonds):8d} !NBOND\n")
        line = []
        for a, b in bonds:
            line.append(f"{a:8d}{b:8d}")
            if len(line) == 4:
                fh.write("".join(line) + "\n")
                line = []
        if line:
            fh.write("".join(line) + "\n")
        fh.write("\n       0 !NIMPHI: impropers\n")
        fh.write("\n       0 !NDON: donors\n")
        fh.write("\n       0 !NACC: acceptors\n")
        fh.write("\n       0 !NNB\n")


def write_pdb(path, atoms, cryst=None):
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        if cryst:
            fh.write("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n" % cryst)
        for i, (name, resname, resid, x, y, z) in enumerate(atoms, 1):
            fh.write("ATOM  %5d %-4s %-3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00\n"
                     % (i, name, resname, resid, x, y, z))
        fh.write("END\n")


def write_xsc(path, step, a=100.0, b=100.0, c=100.0, origin=(50.0, 50.0, 50.0)):
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write("# NAMD extended system configuration restart file\n")
        fh.write("#$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z s_x s_y s_z s_u s_v s_w\n")
        fh.write(f"{step} {a} 0 0 0 {b} 0 0 0 {c} {origin[0]} {origin[1]} {origin[2]} "
                 "0 0 0 0 0 0\n")


class ToolCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = self.tmp.name

    def tearDown(self):
        self.tmp.cleanup()

    def p(self, name):
        return os.path.join(self.dir, name)


class TestDcd(ToolCase):
    def test_roundtrip(self):
        coords = np.array([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0], [-4.0, 5.0, -6.0]], dtype="f4")
        frames = [coords, coords + 1.0, coords + 2.0]
        natoms = write_dcd(self.p("t.dcd"), frames, istart=500, nsavc=250)
        d = namdlib.Dcd(self.p("t.dcd"))
        try:
            self.assertEqual(d.natoms, natoms)
            self.assertEqual(d.frames, 3)
            self.assertEqual(d.istart, 500)
            self.assertEqual(d.nsavc, 250)
            for i in range(3):
                got, cell = d.frame(i)
                np.testing.assert_allclose(got, frames[i], atol=1e-4)
                np.testing.assert_allclose(np.diag(cell), [100.0, 100.0, 100.0])
        finally:
            d.close()

    def test_partial_trailing_frame_is_reported(self):
        coords = np.zeros((4, 3), dtype="f4")
        write_dcd(self.p("t.dcd"), [coords, coords])
        with open(self.p("t.dcd"), "ab") as fh:
            fh.write(b"\0" * 100)
        d = namdlib.Dcd(self.p("t.dcd"))
        try:
            self.assertEqual(d.frames, 2)
            self.assertGreater(d.partial_bytes, 0)
        finally:
            d.close()

    def test_frame_index_subset(self):
        coords = np.arange(30, dtype="f4").reshape(10, 3)
        write_dcd(self.p("t.dcd"), [coords])
        d = namdlib.Dcd(self.p("t.dcd"))
        try:
            got, _ = d.frame(0, indices=np.array([0, 5, 9]))
            np.testing.assert_allclose(got, coords[[0, 5, 9]], atol=1e-4)
        finally:
            d.close()

    def test_frame_bytes_matches_written_size(self):
        coords = np.zeros((7, 3), dtype="f4")
        write_dcd(self.p("t.dcd"), [coords, coords])
        d = namdlib.Dcd(self.p("t.dcd"))
        try:
            self.assertEqual(d.frame_bytes, namdlib.frame_bytes_for(7))
        finally:
            d.close()


class TestOtherReaders(ToolCase):
    def test_psf(self):
        atoms = [(1, "LIG", "C1", "CG2R61", -0.1, 12.011),
                 (1, "LIG", "C2", "CG2R61", 0.1, 12.011),
                 (2, "LIG", "C1", "CG2R61", -0.1, 12.011)]
        write_psf(self.p("a.psf"), atoms, bonds=[(0, 1)])
        n, parsed = namdlib.read_psf(self.p("a.psf"))
        self.assertEqual(n, 3)
        self.assertEqual(parsed[0]["name"], "C1")
        self.assertEqual(parsed[2]["resid"], 2)
        self.assertAlmostEqual(parsed[1]["charge"], 0.1)

    def test_pdb_box(self):
        write_pdb(self.p("a.pdb"), [("C1", "LIG", 1, 1.0, 2.0, 3.0)],
                  cryst=(10.0, 20.0, 30.0, 90.0, 90.0, 90.0))
        box = namdlib.read_pdb_box(self.p("a.pdb"))
        self.assertEqual(box[:3], (10.0, 20.0, 30.0))

    def test_xsc(self):
        write_xsc(self.p("a.xsc"), step=123456, a=11.0, b=12.0, c=13.0)
        xs = namdlib.read_xsc(self.p("a.xsc"))
        self.assertEqual(xs["step"], 123456)
        np.testing.assert_allclose(np.diag(xs["cell"]), [11.0, 12.0, 13.0])


class TestGeometry(ToolCase):
    def test_minimum_image_wraps_into_half_box(self):
        cell = np.diag([10.0, 10.0, 10.0])
        d = namdlib.minimum_image(np.array([9.0, -9.0, 0.0]), cell)
        np.testing.assert_allclose(d, [-1.0, 1.0, 0.0])

    def test_unwrap_molecule_stays_in_global_frame(self):
        """Regression: unwrap used to re-centre each molecule on its own first atom,
        which destroyed all inter-molecular distances."""
        cell = np.diag([100.0, 100.0, 100.0])
        # a molecule straddling the boundary: atoms at x = 99 and x = 1
        mol = np.array([[99.0, 50.0, 50.0], [1.0, 50.0, 50.0]])
        out = namdlib.unwrap_molecule(mol, cell)
        # internally it must be whole: the two atoms 2 A apart, not 98
        self.assertAlmostEqual(float(np.linalg.norm(out[1] - out[0])), 2.0, places=3)
        # and it must stay near its own position, not be moved to the origin
        self.assertGreater(float(np.mean(np.abs(out[:, 0]))), 50.0)

    def test_unwrap_keeps_molecules_apart(self):
        cell = np.diag([100.0, 100.0, 100.0])
        a = namdlib.unwrap_molecule(np.array([[10.0, 10.0, 10.0], [11.0, 10.0, 10.0]]), cell)
        b = namdlib.unwrap_molecule(np.array([[70.0, 70.0, 70.0], [71.0, 70.0, 70.0]]), cell)
        dist = float(np.linalg.norm(a.mean(axis=0) - b.mean(axis=0)))
        self.assertGreater(dist, 30.0)


class TestHelpers(unittest.TestCase):
    def test_frame_bytes_for(self):
        self.assertEqual(namdlib.frame_bytes_for(10), 56 + 3 * (8 + 40))

    def test_trajectory_bytes_scales_with_ns(self):
        one = namdlib.trajectory_bytes(1000, 10.0, 10.0)
        two = namdlib.trajectory_bytes(1000, 20.0, 10.0)
        self.assertAlmostEqual(two / one, 2.0, places=6)

    def test_human_bytes(self):
        self.assertIn("GB", namdlib.human_bytes(5 * 1024 ** 3))


if __name__ == "__main__":
    unittest.main(verbosity=2)