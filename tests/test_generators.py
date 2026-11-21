"""Input generators: box sizing, Packmol/VMD scripts, NAMD .conf rendering."""
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
import preflight  # noqa: E402
import scaffold  # noqa: E402


def sample_case(**over):
    case = {
        "name": "demo",
        "psf": "ionized.psf",
        "pdb": "ionized.pdb",
        "param_dir": "toppar",
        "ligand_param": "lig.prm",
        "ligand_resname": "LIG",
        "forcefield_params": ["par_all36_cgenff.prm", "toppar_water_ions.str"],
        "temperature_K": 300.0,
        "timestep_fs": 2.0,
        "dcdfreq_ps": 10.0,
        "restartfreq_ps": 10.0,
        "energyfreq_ps": 20.0,
        "cutoff": 12.0,
        "switchdist": 10.0,
        "pairlistdist": 14.0,
        "scaling_14": 1.0,
        "stages": {"min": {"minimize_steps": 10000}, "nvt": {"ns": 0.2},
                   "npt": {"ns": 1.0}, "prod": {"ns": 100.0}},
        "box": [[118.0, 0.0, 0.0], [0.0, 118.0, 0.0], [0.0, 0.0, 118.0]],
        "box_origin": [59.0, 59.0, 59.0],
        "natoms": 100000,
        "analysis": {},
    }
    case.update(over)
    return case


class TestBoxMath(unittest.TestCase):
    def test_box_for_concentration(self):
        side, vol = scaffold.box_for_concentration(50, 50.0)
        # 50 molecules at 50 mM -> 1.66e6 A^3 -> 118.4 A
        self.assertAlmostEqual(vol / 1e6, 1.661, places=2)
        self.assertAlmostEqual(side, 118.42, places=1)

    def test_concentration_roundtrip(self):
        """Building a box for c mM and computing c back must return c."""
        for nmol, conc in ((50, 50.0), (100, 25.0), (20, 80.0)):
            side, vol = scaffold.box_for_concentration(nmol, conc)
            back = nmol / (vol * 6.02214076e-4) * 1000.0
            self.assertAlmostEqual(back, conc, places=6)

    def test_box_scales_as_cube_root(self):
        s1, _ = scaffold.box_for_concentration(50, 50.0)
        s8, _ = scaffold.box_for_concentration(400, 50.0)
        self.assertAlmostEqual(s8 / s1, 2.0, places=4)


class TestPackmolScript(unittest.TestCase):
    def test_packing_box_is_smaller_than_final(self):
        text, out = scaffold.render_packmol({"name": "demo", "ligand_pdb": "mol.pdb"},
                                            50, 118.0, pad=8.0)
        self.assertIn("inside box 0. 0. 0. 102.000", text)
        self.assertEqual(out, "demo_com.pdb")

    def test_seed_is_fixed(self):
        text, _ = scaffold.render_packmol({"name": "demo", "ligand_pdb": "mol.pdb"},
                                          50, 118.0, seed=4242)
        self.assertIn("seed 4242", text)


class TestSolvateScript(unittest.TestCase):
    def test_pins_final_box(self):
        text = scaffold.render_solvate({"name": "demo"}, 118.4)
        self.assertIn("-minmax", text)
        self.assertIn("118.400", text)
        self.assertNotIn("-t 8", text)


class TestConfRendering(unittest.TestCase):
    def setUp(self):
        self.case = sample_case()
        self.rendered = scaffold.render_all(self.case)

    def test_all_stages_present(self):
        self.assertEqual(set(self.rendered),
                         {"01_min.conf", "02_nvt.conf", "03_equil_npt.conf", "04_prod.conf"})

    def test_uses_current_keyword_not_deprecated(self):
        for name, text in self.rendered.items():
            self.assertNotIn("1-4scaling", text, f"{name} uses the deprecated keyword")
            self.assertIn("oneFourScaling", text, f"{name} lacks oneFourScaling")

    def test_scaling_value_honoured(self):
        case = sample_case(scaling_14=1.0)
        self.assertIn("oneFourScaling     1.0", scaffold.render("prod", case))

    def test_production_numbers_consistent(self):
        scalars, _ = preflight.parse_conf_dict(
            scaffold.render("prod", self.case)) if hasattr(preflight, "parse_conf_dict") else (None, None)
        if scalars is None:  # no helper; parse from text instead
            import tempfile
            with tempfile.NamedTemporaryFile("w", suffix=".conf", delete=False,
                                             encoding="utf-8") as fh:
                fh.write(scaffold.render("prod", self.case))
                path = fh.name
            try:
                scalars, _ = preflight.parse_conf(path)
            finally:
                os.unlink(path)
        dt = self.case["timestep_fs"] / 1000.0
        self.assertEqual(int(scalars["run"]), int(round(100.0 * 1000.0 / dt)))
        self.assertEqual(int(scalars["dcdfreq"]), int(round(10.0 / dt)))

    def test_min_stage_defines_the_box(self):
        text = scaffold.render("min", self.case)
        self.assertIn("cellBasisVector1", text)
        self.assertIn("minimize           10000", text)

    def test_production_takes_cell_from_xsc(self):
        text = scaffold.render("prod", self.case)
        self.assertIn("extendedSystem", text)
        active = [l for l in text.splitlines()
                  if l.strip().startswith("cellBasisVector")]
        self.assertEqual(active, [], "cellBasisVector must be commented out for a continuation")

    def test_rendered_roundtrips_through_parser(self):
        for name, text in self.rendered.items():
            import tempfile
            with tempfile.NamedTemporaryFile("w", suffix=".conf", delete=False,
                                             encoding="utf-8") as fh:
                fh.write(text)
                path = fh.name
            try:
                scalars, _ = preflight.parse_conf(path)
            finally:
                os.unlink(path)
            self.assertIn("timestep", scalars, name)
            if name != "01_min.conf":
                self.assertIn("run", scalars, name)


if __name__ == "__main__":
    unittest.main(verbosity=2)