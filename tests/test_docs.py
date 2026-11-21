"""Repository hygiene: encoding, syntax, no personal data, CLI entry points."""
import os
import re
import subprocess
import sys
import unittest


def _read(path, encoding="utf-8", errors="replace"):
    with open(path, encoding=encoding, errors=errors) as fh:
        return fh.read()

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
SCRIPTS = os.path.join(ROOT, "scripts")

# Anything that must not appear in a public repository.
FORBIDDEN = re.compile(
    r"C:\\\\Users\\\\|F:\\\\work1|D:\\\\software"      # personal paths
    r"|NpCyJ|npcyj|Npcy"                                # unpublished molecule names
    r"|103\.1|45\.6|1\.58"                              # specific unpublished numbers
)


def walk(exts=None):
    for dirpath, dirnames, files in os.walk(ROOT):
        dirnames[:] = [d for d in dirnames if d not in (".git", "__pycache__", "tests")]
        for f in files:
            if exts and not f.endswith(exts):
                continue
            yield os.path.join(dirpath, f)


class TestEncoding(unittest.TestCase):
    def test_text_files_are_utf8_without_bom(self):
        bad = []
        for p in walk((".md", ".py", ".txt")):
            with open(p, "rb") as fh:
                raw = fh.read()
            if raw.startswith(b"\xef\xbb\xbf"):
                bad.append(f"{os.path.relpath(p, ROOT)}: BOM at start")
                continue
            try:
                txt = raw.decode("utf-8")
            except UnicodeDecodeError:
                bad.append(f"{os.path.relpath(p, ROOT)}: not UTF-8")
                continue
            if "\ufeff" in txt:
                bad.append(f"{os.path.relpath(p, ROOT)}: BOM inside")
        self.assertEqual(bad, [])

    def test_license_is_the_verbatim_gpl(self):
        text = _read(os.path.join(ROOT, "LICENSE"))
        self.assertIn("GNU GENERAL PUBLIC LICENSE", text)
        self.assertIn("Version 3, 29 June 2007", text)
        self.assertIn("END OF TERMS AND CONDITIONS", text)
        # the license text must not be modified
        self.assertNotIn("Copyright (C) 2026", text)


class TestNoPersonalData(unittest.TestCase):
    def test_no_forbidden_strings(self):
        hits = []
        for p in walk((".md", ".py", ".txt")):
            for i, line in enumerate(open(p, encoding="utf-8", errors="replace"), 1):
                if FORBIDDEN.search(line):
                    hits.append(f"{os.path.relpath(p, ROOT)}:{i}: {line.strip()[:70]}")
        self.assertEqual(hits, [])


class TestSyntax(unittest.TestCase):
    def test_every_script_compiles(self):
        for f in sorted(os.listdir(SCRIPTS)):
            if not f.endswith(".py"):
                continue
            path = os.path.join(SCRIPTS, f)
            with open(path, encoding="utf-8") as fh:
                src = fh.read()
            try:
                compile(src, path, "exec")
            except SyntaxError as exc:
                self.fail(f"{f} does not compile: {exc}")

    def test_every_script_has_a_docstring(self):
        for f in sorted(os.listdir(SCRIPTS)):
            if not f.endswith(".py"):
                continue
            src = _read(os.path.join(SCRIPTS, f))
            # module docstring appears after the licence header block
            body = "\n".join(l for l in src.splitlines() if not l.startswith("#"))
            self.assertTrue(body.lstrip().startswith(('"""', "'''")),
                            f"{f} has no module docstring")


class TestEntryPoints(unittest.TestCase):
    """Each tool must answer --help with exit 0, so an install problem surfaces early."""

    TOOLS = ["pipeline", "scaffold", "inspect_molecule", "preflight", "restart",
             "inspect_traj", "aggregate", "analyze", "compare_systems", "study_report"]

    def test_help_exits_zero(self):
        for tool in self.TOOLS:
            path = os.path.join(SCRIPTS, f"{tool}.py")
            r = subprocess.run([sys.executable, path, "--help"],
                               capture_output=True, text=True, timeout=120)
            self.assertEqual(r.returncode, 0, f"{tool} --help failed: {r.stderr[:300]}")

    def test_modules_import(self):
        sys.path.insert(0, SCRIPTS)
        for tool in self.TOOLS:
            __import__(tool)


class TestDocs(unittest.TestCase):
    def test_referenced_docs_exist(self):
        for md in ("README.md", "SKILL.md"):
            text = _read(os.path.join(ROOT, md))
            for target in re.findall(r"\]\((docs/[^)]+\.md)\)", text):
                self.assertTrue(os.path.exists(os.path.join(ROOT, target)),
                                f"{md} links to missing {target}")

    def test_license_mentioned(self):
        text = _read(os.path.join(ROOT, "README.md"))
        self.assertIn("GPL-3.0", text)
        self.assertIn("Copyright (C) 2026", text)


if __name__ == "__main__":
    unittest.main(verbosity=2)