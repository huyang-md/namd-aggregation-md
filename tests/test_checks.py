"""Checks: conf parsing, drift detection, firsttimestep, comparability."""
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
import compare_systems  # noqa: E402
import preflight  # noqa: E402
import scaffold  # noqa: E402


def write(path, text):
    with open(path, "w", encoding="utf-8-sig", newline="\n") as fh:
        fh.write(text)


class TestParseConf(unittest.TestCase):
    def test_comments_and_inline_semicolons(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "a.conf")
            write(p, "timestep 2.0   ;# inline comment\n"
                     "# full line comment\n"
                     "cellBasisVector1   10.0 0.0 0.0\n")
            scalars, repeated = preflight.parse_conf(p)
            self.assertEqual(scalars["timestep"], "2.0")
            self.assertEqual(scalars["cellBasisVector1"], "0.0")
            self.assertEqual(repeated["cellBasisVector1"], ["10.0", "0.0", "0.0"])

    def test_deprecated_scaling_is_aliased(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "a.conf")
            write(p, "1-4scaling 1.07\n")
            scalars, _ = preflight.parse_conf(p)
            self.assertEqual(scalars["oneFourScaling"], "1.07")

    def test_current_scaling_is_left_alone(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "a.conf")
            write(p, "oneFourScaling 1.0\n")
            scalars, _ = preflight.parse_conf(p)
            self.assertEqual(scalars["oneFourScaling"], "1.0")


class TestConsistency(unittest.TestCase):
    def test_flags_differing_scaling(self):
        confs = {
            "a.conf": ({"oneFourScaling": "1.0", "timestep": "2.0"}, {}),
            "b.conf": ({"oneFourScaling": "1.07", "timestep": "2.0"}, {}),
        }
        issues = []
        preflight.check_consistency(".", confs, issues)
        self.assertTrue(any(i["level"] == "FAIL" and "oneFourScaling" in i["message"]
                            for i in issues))

    def test_accepts_identical_settings(self):
        confs = {n: ({"oneFourScaling": "1.0", "timestep": "2.0"}, {})
                 for n in ("a.conf", "b.conf")}
        issues = []
        preflight.check_consistency(".", confs, issues)
        self.assertFalse(any(i["level"] == "FAIL" for i in issues))


class TestDrift(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = self.tmp.name
        self.case = {
            "name": "demo", "psf": "ionized.psf", "pdb": "ionized.pdb",
            "param_dir": "toppar", "ligand_param": "", "ligand_resname": "LIG",
            "forcefield_params": [], "temperature_K": 300.0, "timestep_fs": 2.0,
            "dcdfreq_ps": 10.0, "restartfreq_ps": 10.0, "energyfreq_ps": 20.0,
            "cutoff": 12.0, "switchdist": 10.0, "pairlistdist": 14.0, "scaling_14": 1.0,
            "stages": {"min": {"minimize_steps": 10000}, "nvt": {"ns": 0.2},
                       "npt": {"ns": 1.0}, "prod": {"ns": 100.0}},
            "box": [[100.0, 0, 0], [0, 100.0, 0], [0, 0, 100.0]],
            "box_origin": [50.0, 50.0, 50.0], "analysis": {},
        }

    def tearDown(self):
        self.tmp.cleanup()

    def test_generated_file_in_sync(self):
        for name, text in scaffold.render_all(self.case).items():
            write(os.path.join(self.dir, name), text)
        issues = []
        preflight.check_drift(self.dir, self.case, {}, issues)
        self.assertFalse(any(i["level"] == "FAIL" for i in issues))

    def test_generated_file_edited_is_flagged(self):
        rendered = scaffold.render_all(self.case)
        text = rendered["04_prod.conf"].replace("oneFourScaling     1.0",
                                                "oneFourScaling     1.07")
        write(os.path.join(self.dir, "04_prod.conf"), text)
        for name, body in rendered.items():
            if name != "04_prod.conf":
                write(os.path.join(self.dir, name), body)
        issues = []
        preflight.check_drift(self.dir, self.case, {}, issues)
        self.assertTrue(any(i["level"] == "FAIL" and i["code"] == "drift" for i in issues))

    def test_hand_written_conf_is_not_drift(self):
        """Adopting a case.json into an existing project must not raise drift FAILs."""
        for name, text in scaffold.render_all(self.case).items():
            body = "# hand written by a human\n" + "\n".join(
                l for l in text.splitlines() if not l.startswith("# GENERATED"))
            write(os.path.join(self.dir, name), body)
        issues = []
        preflight.check_drift(self.dir, self.case, {}, issues)
        self.assertFalse(any(i["level"] == "FAIL" for i in issues))
        self.assertTrue(any(i["code"] == "unmanaged-conf" for i in issues))


class TestFirstTimestep(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = self.tmp.name

    def tearDown(self):
        self.tmp.cleanup()

    def _xsc(self, step):
        with open(os.path.join(self.dir, "03_equil_npt.xsc"), "w", encoding="ascii") as fh:
            fh.write("# header\n")
            fh.write(f"{step} 100 0 0 0 100 0 0 0 100 50 50 50 0 0 0 0 0 0\n")

    def test_mismatch_is_fail(self):
        self._xsc(500000)
        issues = []
        preflight.check_static("04_prod.conf",
                               {"extendedSystem": "03_equil_npt.xsc", "firsttimestep": "123"},
                               self.dir, issues)
        self.assertTrue(any(i["level"] == "FAIL" and i["code"] == "firsttimestep-mismatch"
                            for i in issues))

    def test_zero_reset_is_only_a_warning(self):
        self._xsc(500000)
        issues = []
        preflight.check_static("04_prod.conf",
                               {"extendedSystem": "03_equil_npt.xsc", "firsttimestep": "0"},
                               self.dir, issues)
        self.assertFalse(any(i["level"] == "FAIL" for i in issues))
        self.assertTrue(any(i["code"] == "firsttimestep-reset" for i in issues))

    def test_temperature_with_binvelocities_is_fail(self):
        issues = []
        preflight.check_static("04_prod.conf",
                               {"binVelocities": "x.vel", "temperature": "300"},
                               self.dir, issues)
        self.assertTrue(any(i["code"] == "temperature+binVelocities" for i in issues))


class TestComparability(unittest.TestCase):
    def _row(self, name, conc, spacing=10.0, ns=100.0, scaling=1.0):
        return {"name": name, "conc_mM": conc, "spacing_ps": [spacing],
                "prod_ns": ns, "scaling_14": scaling, "molecules": 50}

    def test_equal_systems_pass(self):
        issues = compare_systems.compare([self._row("a", 50.0), self._row("b", 50.0)])
        self.assertFalse(any(level == "FAIL" for level, _k, _m in issues))

    def test_concentration_spread_blocks(self):
        issues = compare_systems.compare([self._row("a", 54.0), self._row("b", 34.0)])
        self.assertTrue(any(k == "concentration" and level == "FAIL" for level, k, _m in issues))

    def test_spacing_spread_blocks(self):
        issues = compare_systems.compare([self._row("a", 50.0, spacing=10.0),
                                          self._row("b", 50.0, spacing=100.0)])
        self.assertTrue(any(k == "spacing" and level == "FAIL" for level, k, _m in issues))

    def test_scaling_spread_blocks(self):
        issues = compare_systems.compare([self._row("a", 50.0, scaling=1.0),
                                          self._row("b", 50.0, scaling=1.07)])
        self.assertTrue(any(k == "scaling_14" and level == "FAIL" for level, k, _m in issues))


if __name__ == "__main__":
    unittest.main(verbosity=2)