"""Output artefacts: the VMD colouring script and the optional PNG plots."""
import os
import sys
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)                                   # reuse the synthetic builders
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))

import aggregate  # noqa: E402
from test_core import ToolCase, write_dcd, write_psf  # noqa: E402

try:
    import matplotlib  # noqa: F401
    HAVE_MPL = True
except ImportError:
    HAVE_MPL = False


def _three_molecule_case(d, add_water=False):
    """Two molecules touching at 5 A centre distance, one far away.

    Returns the case directory; run aggregate.main on it and the first two
    molecules must come out as one cluster, the third as a singleton.
    """
    centers = [(50.0, 50.0, 50.0), (55.0, 50.0, 50.0), (20.0, 20.0, 20.0)]
    atoms = []
    for resid, center in enumerate(centers, 1):
        for k in range(8):
            atoms.append((resid, "DYE", f"C{k + 1}", "C", 0.0, 12.010))
    if add_water:
        atoms.append((4, "TIP3", "OH2", "O", 0.0, 16.0))
    offs = np.array([(sx * 1.25, sy * 1.25, sz * 1.25)
                     for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)])
    frames = [np.vstack([np.asarray(c) + offs for c in centers]) for _ in range(2)]
    if add_water:
        frames = [np.vstack([f, [[10.0, 10.0, 10.0]]]) for f in frames]
    write_psf(os.path.join(d, "ionized.psf"), atoms)
    write_dcd(os.path.join(d, "t.dcd"), frames)
    return d


class TestVmdScript(ToolCase):
    def test_clusters_and_singletons(self):
        d = _three_molecule_case(self.dir)
        rc = aggregate.main(["--dir", d, "--traj", "t.dcd", "--vmd-out"])
        self.assertEqual(rc, 0)
        path = os.path.join(d, "analysis", "clusters_vmd.tcl")
        with open(path, encoding="utf-8") as fh:
            tcl = fh.read()
        self.assertIn("mol delrep 0", tcl)
        # the touching pair share one coloured selection
        self.assertIn("resname DYE and resid 1 2", tcl)
        # the far molecule ends up in the grey singleton selection
        self.assertIn("resid 3", tcl)
        # no water in the PSF -> no solvent representation, 2 addreps total
        self.assertNotIn("water", tcl)
        self.assertEqual(tcl.count("mol addrep"), 2)

    def test_solvent_gets_its_own_dim_rep(self):
        d = _three_molecule_case(self.dir, add_water=True)
        rc = aggregate.main(["--dir", d, "--traj", "t.dcd", "--vmd-out"])
        self.assertEqual(rc, 0)
        path = os.path.join(d, "analysis", "clusters_vmd.tcl")
        with open(path, encoding="utf-8") as fh:
            tcl = fh.read()
        self.assertIn("water or ion", tcl)
        # solvent rep + one cluster rep + one singleton rep
        self.assertEqual(tcl.count("mol addrep"), 3)


@unittest.skipUnless(HAVE_MPL, "matplotlib not installed")
class TestPlots(ToolCase):
    def test_all_pngs_are_rendered(self):
        size_counts = {1: 4, 2: 2, 5: 1}
        series = [(0.0, 3, 5, 1.6), (0.2, 3, 5, 1.7), (0.4, 2, 4, 1.8)]
        shape_rows = [(0.0, 5, 6.1, 0.1, 0.1, 0.05, "sphere"),
                      (0.2, 5, 6.2, 0.1, 0.1, 0.05, "sphere"),
                      (0.4, 4, 5.9, 0.1, 0.1, 0.05, "sphere")]
        shell = np.arange(0.0, 10.0, 1.0)
        dens = np.ones(len(shell) - 1)
        n = aggregate.write_plots(self.dir, size_counts, 7, series, shape_rows,
                                  (shell, dens), "time (ns)")
        self.assertEqual(n, 4)
        for name in ("clusters_distribution.png", "clusters_timeseries.png",
                     "shape_rg.png", "density_profile.png"):
            p = os.path.join(self.dir, name)
            self.assertTrue(os.path.exists(p) and os.path.getsize(p) > 1000, name)

    def test_without_density_profile(self):
        n = aggregate.write_plots(self.dir, {1: 1, 3: 1}, 2,
                                  [(0.0, 1, 3, 2.0), (1.0, 1, 3, 2.0)],
                                  [(0.0, 3, 5.0, 0.1, 0.1, 0.05, "sphere"),
                                   (1.0, 3, 5.1, 0.1, 0.1, 0.05, "sphere")],
                                  None, "frame")
        self.assertEqual(n, 3)


if __name__ == "__main__":
    unittest.main(verbosity=2)
