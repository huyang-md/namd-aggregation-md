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

"""Pre-flight checks for a NAMD/CHARMM-GUI case directory.

Catches, before a long run starts:
  * hand-edits that drifted from case.json
  * inconsistent 1-4scaling / cutoff / timestep / dcdfreq across .conf files
  * restart files whose step disagrees with firsttimestep
  * cellBasisVector set on top of extendedSystem
  * temperature set together with binVelocities (NAMD stops with FATAL)
  * a production run that would overwrite an existing trajectory
  * not enough free disk for the requested trajectory

Exit code: 0 = no FAIL, 2 = at least one FAIL.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import namdlib  # noqa: E402
import scaffold  # noqa: E402

FAIL, WARN, INFO = "FAIL", "WARN", "INFO"
KEY_FIELDS = ("oneFourScaling", "cutoff", "timestep", "exclude", "rigidBonds")


def parse_conf(path):
    """Parse a NAMD .conf into (scalars, repeated-keys)."""
    scalars, repeated = {}, defaultdict(list)
    with open(path, "r", encoding="utf-8-sig", errors="replace") as fh:
        for raw in fh:
            line = raw.split("#", 1)[0]
            if ";" in line:
                line = line.split(";", 1)[0]
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            key = parts[0]
            if key == "set" and len(parts) >= 3:
                scalars[f"${parts[1]}"] = parts[2]
                continue
            vals = parts[1:]
            if not vals:
                continue
            repeated[key].extend(vals)
            scalars[key] = vals[-1]
    # "1-4scaling" is deprecated in current NAMD (it warns and asks for
    # oneFourScaling); treat both spellings as the same knob.
    if "oneFourScaling" not in scalars and "1-4scaling" in scalars:
        scalars["oneFourScaling"] = scalars["1-4scaling"]
    return scalars, repeated


def add(issues, level, code, message):
    issues.append(dict(level=level, code=code, message=message))


def check_consistency(directory, confs, issues):
    """Cross-file consistency of physics-relevant keywords."""
    seen = defaultdict(list)
    for name, (scalars, _) in confs.items():
        for field in KEY_FIELDS:
            if field in scalars:
                seen[field].append((name, scalars[field]))
    for field, entries in seen.items():
        values = {v for _, v in entries}
        if len(values) > 1:
            detail = ", ".join(f"{n}={v}" for n, v in sorted(entries))
            add(issues, FAIL, f"inconsistent:{field}",
                f"{field} differs between .conf files -> {detail}")
        else:
            add(issues, INFO, f"consistent:{field}", f"{field} = {entries[0][1]} in all confs")


def check_static(name, scalars, directory, issues, known_outputs=()):
    """Single-file sanity checks."""
    bt = scalars.get("binVelocities")
    if bt and "temperature" in scalars and "reassignHold" not in scalars:
        add(issues, FAIL, "temperature+binVelocities",
            f"{name}: 'temperature' is set together with binVelocities {bt}; "
            "NAMD aborts. Remove the temperature line.")

    es = scalars.get("extendedSystem")
    if es and any(k.startswith("cellBasisVector") for k in scalars):
        add(issues, WARN, "cell+extendedSystem",
            f"{name}: cellBasisVector is set while extendedSystem {es} is used; "
            "the cell should come from the .xsc only.")

    if "1-4scaling" in scalars:
        add(issues, WARN, "deprecated-keyword",
            f"{name}: '1-4scaling' is deprecated in current NAMD; it still works but prints "
            "a warning. Rename it to 'oneFourScaling', which is what NAMD actually honours.")

    for key in ("binCoordinates", "binVelocities", "extendedSystem"):
        f = scalars.get(key)
        if not f:
            continue
        if os.path.exists(os.path.join(directory, f)):
            continue
        stem = os.path.basename(f).split(".")[0]
        if stem in known_outputs:
            add(issues, INFO, f"pending:{key}",
                f"{name}: {f} is the output of stage '{stem}', not produced yet")
        else:
            add(issues, FAIL, f"missing:{key}", f"{name}: {key} points to missing {f}")

    ft = scalars.get("firsttimestep")
    es = scalars.get("extendedSystem")
    if es and ft is not None and os.path.exists(os.path.join(directory, es)):
        try:
            step = namdlib.read_xsc(os.path.join(directory, es))["step"]
        except Exception as exc:  # noqa: BLE001
            add(issues, WARN, "xsc-unreadable", f"{name}: cannot read {es} ({exc})")
            return
        if int(float(ft)) != step:
            if int(float(ft)) == 0:
                add(issues, WARN, "firsttimestep-reset",
                    f"{name}: firsttimestep 0 while {es} is at step {step}; the output "
                    f"numbering restarts at 0, so a later continuation of this stage needs "
                    f"firsttimestep {step}.")
            else:
                add(issues, FAIL, "firsttimestep-mismatch",
                    f"{name}: firsttimestep={int(float(ft))} but {es} is at step {step}. "
                    f"Use firsttimestep {step}.")

    out = scalars.get("outputName")
    if out:
        dcd = os.path.join(directory, out + ".dcd")
        if os.path.exists(dcd) and os.path.getsize(dcd) > 0:
            add(issues, WARN, "overwrite-risk",
                f"{name}: outputName {out} would append to the existing {out}.dcd "
                f"({namdlib.human_bytes(os.path.getsize(dcd))}). "
                "Use a new outputName for a restart.")


def check_drift(directory, case, confs, issues):
    if case is None:
        add(issues, INFO, "no-case",
            "no case.json here. Run 'scaffold.py --dir <dir>' WITHOUT --force to add one: "
            "existing .conf files are left untouched and only case.json is written.")
        return
    rendered = scaffold.render_all(case)

    # A directory is "managed" only if at least one .conf carries the generator header.
    # Otherwise it is a hand-written project that merely adopted a case.json, and the
    # canonical stage names may legitimately differ (min.conf vs 01_min.conf).
    managed = False
    for f in os.listdir(directory):
        if not f.endswith(".conf"):
            continue
        try:
            with open(os.path.join(directory, f), encoding="utf-8-sig", errors="replace") as fh:
                if "GENERATED by namd-aggregation-md" in fh.read(800):
                    managed = True
                    break
        except OSError:
            continue

    for name, text in rendered.items():
        path = os.path.join(directory, name)
        if not os.path.exists(path):
            if managed:
                add(issues, WARN, "missing-conf", f"{name} from case.json is not on disk")
            else:
                add(issues, INFO, "unmanaged-conf",
                    f"no {name}; this directory uses its own stage names, which is fine")
            continue
        disk = open(path, encoding="utf-8-sig", errors="replace").read().replace("\r\n", "\n")
        if disk.strip() == text.strip():
            add(issues, INFO, "in-sync", f"{name} matches case.json")
            continue
        # Only a file this tool generated can have "drifted". A hand-written conf in a
        # directory that merely adopted a case.json is the expected state, not a fault.
        head = "\n".join(disk.splitlines()[:8])
        if "GENERATED by namd-aggregation-md" not in head:
            add(issues, INFO, "unmanaged-conf",
                f"{name} is hand-written, so the drift check does not apply "
                "(case.json still supplies natoms/box/analysis settings)")
            continue
        dl = [l for l in disk.splitlines() if l.strip() and not l.lstrip().startswith("#")]
        rl = [l for l in text.splitlines() if l.strip() and not l.lstrip().startswith("#")]
        diff = _first_diffs(dl, rl)
        add(issues, FAIL, "drift", f"{name} differs from case.json: {diff}")


def _first_diffs(a, b, limit=3):
    import difflib
    out = []
    for line in difflib.unified_diff(a, b, lineterm="", n=0):
        if line.startswith(("---", "+++", "@@")):
            continue
        out.append(line.strip())
        if len(out) >= limit:
            break
    return " | ".join(out) if out else "(whitespace only)"


def check_disk(directory, case, confs, issues):
    if case is not None:
        natoms = case.get("natoms")
        ns = case["stages"]["prod"]["ns"]
        psf = os.path.join(directory, case["psf"])
        if not natoms and os.path.exists(psf):
            natoms, _ = namdlib.read_psf(psf)
        if natoms:
            need = namdlib.trajectory_bytes(natoms, ns, case["dcdfreq_ps"])
            _report_disk(directory, need, f"prod {ns} ns @ {case['dcdfreq_ps']} ps", issues)
        return
    name = scaffold.STAGE_FILES["prod"]
    if name not in confs:
        return
    scalars, _ = confs[name]
    psf = scalars.get("structure")
    if not psf or not os.path.exists(os.path.join(directory, psf)):
        return
    natoms, _ = namdlib.read_psf(os.path.join(directory, psf))
    try:
        steps = float(scalars.get("run", "0"))
        dcd = float(scalars.get("dcdfreq", "0"))
        ts = float(scalars.get("timestep", "2.0"))
    except ValueError:
        return
    if dcd <= 0:
        return
    frames = steps / dcd
    need = int(frames * namdlib.frame_bytes_for(natoms))
    ns = steps * ts / 1e6
    _report_disk(directory, need, f"prod {ns:.1f} ns ({frames:.0f} frames)", issues)


def _report_disk(directory, need, label, issues):
    try:
        free = shutil.disk_usage(directory).free
    except OSError as exc:
        add(issues, WARN, "disk-unknown", f"cannot read free space: {exc}")
        return
    msg = f"{label} needs ~{namdlib.human_bytes(need)}, free on this drive: {namdlib.human_bytes(free)}"
    if need > free * 0.9:
        add(issues, FAIL, "disk", msg + " -> not enough, the run will hit 'No space left on device'")
    elif need > free * 0.5:
        add(issues, WARN, "disk", msg + " -> tight")
    else:
        add(issues, INFO, "disk", msg)


def check_inputs(directory, case, confs, issues):
    if case:
        psf, pdb = case["psf"], case["pdb"]
    else:
        ref = confs.get(scaffold.STAGE_FILES["min"]) or next(iter(confs.values()))
        psf, pdb = ref[0].get("structure"), ref[0].get("coordinates")
    for label, f in (("structure", psf), ("coordinates", pdb)):
        if f and not os.path.exists(os.path.join(directory, f)):
            add(issues, FAIL, "missing-input", f"{label} file {f} is missing")
    if psf and os.path.exists(os.path.join(directory, psf)):
        natoms, atoms = namdlib.read_psf(os.path.join(directory, psf))
        add(issues, INFO, "natoms", f"{psf}: {natoms:,} atoms")
        if case and case.get("natoms") and case["natoms"] != natoms:
            add(issues, FAIL, "natoms-mismatch",
                f"case.json says {case['natoms']:,} atoms, {psf} has {natoms:,}")
        resname = (case or {}).get("ligand_resname")
        if resname:
            nlig = sum(1 for a in atoms if a["resname"] == resname)
            nmol = len({a["resid"] for a in atoms if a["resname"] == resname})
            if nlig:
                add(issues, INFO, "ligand",
                    f"{nmol} x {resname} residues, {nlig:,} atoms "
                    f"({nlig // nmol if nmol else 0} atoms/molecule)")
    if pdb and os.path.exists(os.path.join(directory, pdb)):
        cryst = namdlib.read_pdb_box(os.path.join(directory, pdb))
        if cryst:
            add(issues, INFO, "cryst1", "CRYST1 %g x %g x %g" % cryst[:3])
            ref = confs.get(scaffold.STAGE_FILES["min"])
            rep = ref[1] if ref else {}
            if all(f"cellBasisVector{i}" in rep for i in (1, 2, 3)):
                try:
                    vecs = [[float(x) for x in rep[f"cellBasisVector{i}"][:3]] for i in (1, 2, 3)]
                    conf_len = [(x * x + y * y + z * z) ** 0.5 for x, y, z in vecs]
                except (ValueError, IndexError):
                    conf_len = None
                if conf_len:
                    dmax = max(abs(a - b) for a, b in zip(conf_len, cryst[:3]))
                    lvl = WARN if dmax > 0.05 else INFO
                    add(issues, lvl, "box",
                        f"min.conf box lengths {[round(x, 3) for x in conf_len]} vs "
                        f"CRYST1 {list(cryst[:3])} (max diff {dmax:.4f} A)")


def check_params(directory, confs, issues):
    """Every `parameters` line must resolve, or NAMD stops at startup."""
    for name, (_scalars, repeated) in confs.items():
        listed = repeated.get("parameters", [])
        if not listed:
            continue
        missing = [p for p in listed if not os.path.exists(os.path.join(directory, p))]
        if missing:
            shown = ", ".join(missing[:3]) + (" ..." if len(missing) > 3 else "")
            add(issues, FAIL, "missing-params",
                f"{name}: {len(missing)} of {len(listed)} parameter file(s) not found: {shown}")
        else:
            add(issues, INFO, "params-ok", f"{name}: all {len(listed)} parameter files present")


def run_checks(directory, verbose=False, recursive=False):
    issues = []
    case = None
    cpath = os.path.join(directory, scaffold.CASE_FILE)
    if os.path.exists(cpath):
        try:
            case = scaffold.load_case(directory)
        except Exception as exc:  # noqa: BLE001
            add(issues, FAIL, "case-unreadable", f"case.json: {exc}")

    confs = {}
    if recursive:
        walker = ((r, fs) for r, _d, fs in os.walk(directory))
    else:
        walker = [(directory, [f for f in os.listdir(directory)
                               if os.path.isfile(os.path.join(directory, f))])]
    for root, files in walker:
        for f in files:
            if f.endswith(".conf"):
                p = os.path.join(root, f)
                try:
                    confs[os.path.relpath(p, directory)] = parse_conf(p)
                except Exception as exc:  # noqa: BLE001
                    add(issues, WARN, "conf-unreadable", f"{f}: {exc}")
    if not confs:
        add(issues, FAIL, "no-conf", "no .conf files in this directory")
        return issues

    check_inputs(directory, case, confs, issues)
    check_params(directory, confs, issues)
    check_consistency(directory, confs, issues)
    known = set(scaffold.STAGE_OUT.values())
    for _scalars, _rep in confs.values():
        if _scalars.get("outputName"):
            known.add(_scalars["outputName"])
    for name, (scalars, _) in confs.items():
        if name == "case.json":
            continue
        check_static(name, scalars, directory, issues, known)
    check_drift(directory, case, confs, issues)
    check_disk(directory, case, confs, issues)
    return issues


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dir", required=True)
    p.add_argument("--json", action="store_true")
    p.add_argument("--recursive", action="store_true",
                   help="also inspect .conf files in subdirectories")
    p.add_argument("--verbose", action="store_true", help="also print INFO lines")
    a = p.parse_args(argv)

    directory = os.path.abspath(a.dir)
    if not os.path.isdir(directory):
        sys.exit(f"not a directory: {directory}")
    issues = run_checks(directory, a.verbose, a.recursive)

    if a.json:
        print(json.dumps(dict(directory=directory, issues=issues), indent=2, ensure_ascii=False))
    else:
        order = {FAIL: 0, WARN: 1, INFO: 2}
        for issue in sorted(issues, key=lambda i: order[i["level"]]):
            if issue["level"] == INFO and not a.verbose:
                continue
            tag = {"FAIL": "[FAIL]", "WARN": "[warn]", "INFO": "[info]"}[issue["level"]]
            print(f"{tag:7} {issue['message']}")
        n_fail = sum(1 for i in issues if i["level"] == FAIL)
        n_warn = sum(1 for i in issues if i["level"] == WARN)
        print(f"\n{directory}: {n_fail} FAIL, {n_warn} WARN, "
              f"{sum(1 for i in issues if i['level'] == INFO)} info")
        if n_fail == 0:
            print("ready to run.")

    return 2 if any(i["level"] == FAIL for i in issues) else 0


if __name__ == "__main__":
    sys.exit(main())
