"""Analysis maths: shape descriptors, periodic-boundary handling, clustering."""
import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
import aggregate  # noqa: E402
import analyze  # noqa: E402


class TestGyration(unittest.TestCase):
    def test_single_atom_is_degenerate(self):
        rg, asp, acyl, kappa, shape = aggregate.gyration(np.zeros((1, 3)))
        self.assertEqual(shape, "single")
        self.assertEqual(rg, 0.0)

    def test_uniform_sphere_is_isotropic(self):
        rng = np.random.default_rng(0)
        # points on a sphere surface are isotropic by symmetry
        v = rng.normal(size=(2000, 3))
        v /= np.linalg.norm(v, axis=1)[:, None]
        rg, asp, acyl, kappa, shape = aggregate.gyration(v * 10.0)
        self.assertLess(kappa, 0.05)
        self.assertEqual(shape, "sphere")

    def test_linear_chain_is_a_rod(self):
        t = np.linspace(-20, 20, 200)
        pts = np.column_stack([t, np.zeros_like(t), np.zeros_like(t)])
        rg, asp, acyl, kappa, shape = aggregate.gyration(pts)
        self.assertGreater(kappa, 0.5)
        self.assertEqual(shape, "rod")
        # compare against the definition rather than the continuum limit
        c = pts.mean(axis=0)
        self.assertAlmostEqual(rg, float(np.sqrt(((pts - c) ** 2).sum(axis=1).mean())), places=6)

    def test_flat_slab_is_discriminated(self):
        rng = np.random.default_rng(1)
        xy = rng.uniform(-10, 10, size=(3000, 2))
        pts = np.column_stack([xy, np.zeros(len(xy))])
        rg, asp, acyl, kappa, shape = aggregate.gyration(pts)
        self.assertIn(shape, ("slab", "irregular"))
        self.assertGreater(kappa, 0.1)

    def test_rg_matches_definition(self):
        pts = np.array([[0.0, 0, 0], [2.0, 0, 0], [0, 2.0, 0], [0, 0, 2.0]])
        c = pts.mean(axis=0)
        expected = float(np.sqrt(((pts - c) ** 2).sum(axis=1).mean()))
        rg, *_ = aggregate.gyration(pts)
        self.assertAlmostEqual(rg, expected, places=6)


class TestPeriodic(unittest.TestCase):
    """Regression: a cluster straddling the boundary must not get a huge Rg."""

    def test_gather_cluster_rejoins_across_boundary(self):
        cell = np.diag([100.0, 100.0, 100.0])
        # two molecules, one at x=2 and one at x=98: 4 A apart through the boundary
        m0 = np.array([[2.0, 50.0, 50.0], [3.0, 50.0, 50.0]])
        m1 = np.array([[98.0, 50.0, 50.0], [99.0, 50.0, 50.0]])
        pts = aggregate.gather_cluster([m0, m1], [0, 1], cell)
        span = pts.max(axis=0) - pts.min(axis=0)
        self.assertLess(float(span[0]), 10.0, "molecules were not brought together")
        rg, *_ = aggregate.gyration(pts)
        self.assertLess(rg, 10.0)

    def test_gather_cluster_keeps_separate_clusters_apart(self):
        cell = np.diag([100.0, 100.0, 100.0])
        m0 = np.array([[10.0, 10.0, 10.0]])
        m1 = np.array([[60.0, 60.0, 60.0]])
        pts = aggregate.gather_cluster([m0, m1], [0, 1], cell)
        self.assertGreater(float(np.linalg.norm(pts[1] - pts[0])), 30.0)


class TestClustering(unittest.TestCase):
    def test_two_groups_are_separated(self):
        cell = np.diag([100.0, 100.0, 100.0])
        # molecules 0,1 close together; 2,3 close together; groups far apart
        mols = [
            np.array([[10.0, 10.0, 10.0]]),
            np.array([[13.0, 10.0, 10.0]]),
            np.array([[60.0, 60.0, 60.0]]),
            np.array([[63.0, 60.0, 60.0]]),
        ]
        labels = aggregate.cluster_labels(mols, cell, cutoff=4.0, mode="contact")
        self.assertEqual(labels[0], labels[1])
        self.assertEqual(labels[2], labels[3])
        self.assertNotEqual(labels[0], labels[2])

    def test_centroid_mode_respects_cutoff(self):
        cell = np.diag([100.0, 100.0, 100.0])
        mols = [np.array([[10.0, 10.0, 10.0]]), np.array([[16.0, 10.0, 10.0]])]
        near = aggregate.cluster_labels(mols, cell, cutoff=8.0, mode="centroid")
        far = aggregate.cluster_labels(mols, cell, cutoff=3.0, mode="centroid")
        self.assertEqual(near[0], near[1])
        self.assertNotEqual(far[0], far[1])

    def test_groups_sorted_by_size(self):
        cell = np.diag([200.0, 200.0, 200.0])
        mols = [np.array([[10.0, 10, 10]]), np.array([[12.0, 10, 10]]),
                np.array([[150.0, 150, 150]])]
        labels = aggregate.cluster_labels(mols, cell, cutoff=4.0, mode="contact")
        groups = aggregate.group_by_label(labels)
        self.assertEqual(len(groups), 2)
        self.assertEqual(len(groups[0]), 2)


class TestAnalyzeHelpers(unittest.TestCase):
    def test_safe_name_strips_path_separators(self):
        self.assertEqual(analyze.safe_name("ring:ring1"), "ring-ring1")
        self.assertEqual(analyze.safe_name("C5|O3"), "C5-O3")

    def test_plane_normal_of_a_known_plane(self):
        # points in the z = 5 plane -> normal along z
        pts = np.array([[0.0, 0, 5], [1.0, 0, 5], [0, 1.0, 5], [1.0, 1.0, 5]])
        c, n, rms = analyze.namdlib.plane(pts) if hasattr(analyze, "namdlib") else (None, None, None)
        if n is None:
            import namdlib
            c, n, rms = namdlib.plane(pts)
        self.assertAlmostEqual(abs(float(n[2])), 1.0, places=6)
        self.assertLess(rms, 1e-6)


if __name__ == "__main__":
    unittest.main(verbosity=2)