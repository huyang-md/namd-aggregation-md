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

"""Generic self-assembly observables — works for any solute, not just chromophores.

Answers the questions every self-assembly study asks, whatever the chemistry:
how big are the aggregates, what shape are they, what are they made of, and how
are they organised inside?

    python aggregate.py --dir <case> --traj 04_prod.dcd --last-ns 50
    python aggregate.py --dir <case> --traj nowater.dcd --frame-ps 100 --mode contact

Complements analyze.py: that one handles stacking geometry (planes, slip, J/H);
this one handles aggregation statistics and structure, which apply to micelles,
vesicles, fibrils, co-assemblies and inorganic clusters alike.

Outputs (CSV, Origin/Grace friendly):
    clusters.csv              size distribution over the sampled frames
    clusters_timeseries.csv   per-frame cluster count, largest size, mean size
    shape.csv                 per-frame size, Rg and shape descriptors of the largest cluster
    composition.csv           molecule-type make-up of each cluster (multi-component only)
    density_profile.csv       radial density of solute around the largest cluster's centre
    summary.txt               human summary
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import Counter, defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import namdlib  # noqa: E402
import scaffold  # noqa: E402

WATER = {"TIP3", "TIP2", "WAT", "HOH", "SOL", "TIPS", "SPC", "TIP4"}
IONS = {"SOD", "CLA", "POT", "CAL", "MG", "ZN", "NA", "CL", "K", "CA", "CES", "BRO", "IOD"}


# --------------------------------------------------------------------------- #
# molecule bookkeeping
# --------------------------------------------------------------------------- #
def solute_molecules(atoms):
    """Group solute atoms into molecules: list of (resname, [atom dicts])."""
    mols = defaultdict(list)
    for a in atoms:
        rn = a["resname"]
        if rn in WATER or rn in IONS:
            continue
        mols[(rn, a["resid"])].append(a)
    ordered = [((rn, rid), mols[(rn, rid)]) for rn, rid in sorted(mols)]
    return ordered


def index_map(dcd, psf_natoms, molecules):
    """Return (indices into the trajectory, atoms-per-molecule list, description)."""
    per = [len(ats) for _key, ats in molecules]
    total = sum(per)
    if dcd.natoms == psf_natoms:
        idx = np.concatenate([np.array([a["serial"] - 1 for a in ats]) for _k, ats in molecules])
        return idx, per, "full system (indices from PSF)"
    if dcd.natoms == total:
        return np.arange(total), per, "solute-only trajectory"
    raise SystemExit(
        f"trajectory has {dcd.natoms:,} atoms; PSF has {psf_natoms:,} and the solute "
        f"totals {total:,}. Neither matches — pass the matching PSF.")


# --------------------------------------------------------------------------- #
# clustering
# --------------------------------------------------------------------------- #
def cluster_labels(coords_by_mol, cell, cutoff, mode):
    """Union-find over molecules; returns a label array (one per molecule)."""
    n = len(coords_by_mol)
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[rj] = ri

    cents = np.array([c.mean(axis=0) for c in coords_by_mol])
    d = cents[:, None, :] - cents[None, :, :]
    d = namdlib.minimum_image(d, cell)
    cdist = np.linalg.norm(d, axis=2)
    rows, cols = np.where(np.triu(cdist < cutoff + (12.0 if mode == "contact" else 0.0), k=1))
    for i, j in zip(rows, cols):
        if mode == "centroid":
            if cdist[i, j] <= cutoff:
                union(i, j)
            continue
        delta = coords_by_mol[i][:, None, :] - coords_by_mol[j][None, :, :]
        delta = namdlib.minimum_image(delta, cell)
        if float(np.sqrt((delta ** 2).sum(axis=2)).min()) <= cutoff:
            union(i, j)
    return np.array([find(i) for i in range(n)])


def group_by_label(labels):
    groups = defaultdict(list)
    for i, lab in enumerate(labels):
        groups[lab].append(i)
    return sorted(groups.values(), key=len, reverse=True)


# --------------------------------------------------------------------------- #
# shape
# --------------------------------------------------------------------------- #
def gather_cluster(mol_xyz, members, cell):
    """Collect a cluster's atoms in a mutually consistent frame.

    Each molecule is already whole internally, but two molecules in the same cluster can
    still sit on opposite sides of the periodic boundary. Shift every molecule by the
    minimum-image difference of its centroid from the first member's centroid.
    """
    ref = mol_xyz[members[0]].mean(axis=0)
    out = []
    for m in members:
        c = mol_xyz[m].mean(axis=0)
        shift = ref + namdlib.minimum_image(c - ref, cell) - c
        out.append(mol_xyz[m] + shift)
    return np.concatenate(out)


def gyration(atoms_xyz):
    """Return (Rg, asphericity, acylindricity, relative anisotropy, shape class)."""
    c = atoms_xyz.mean(axis=0)
    r = atoms_xyz - c
    # gyration tensor (mass-weighted would need masses; atom count is a fine proxy)
    g = (r.T @ r) / len(atoms_xyz)
    ev = np.sort(np.linalg.eigvalsh(g))          # ascending
    ev = np.clip(ev, 0.0, None)
    rg = float(np.sqrt(ev.sum()))
    if rg <= 0:
        return 0.0, 0.0, 0.0, 0.0, "single"
    l1, l2, l3 = ev
    asphericity = float(l3 - 0.5 * (l1 + l2)) / (rg ** 2)
    acylindricity = float(l2 - l1) / (rg ** 2)
    trace2 = float(l1 + l2 + l3) ** 2
    kappa2 = 1.0 - 3.0 * float(l1 * l2 + l2 * l3 + l3 * l1) / max(trace2, 1e-12)
    # crude classification from the eigenvalue ratios
    ratios = ev / max(l3, 1e-12)
    if kappa2 < 0.05:
        shape = "sphere"
    elif ratios[0] < 0.15 and ratios[1] < 0.35:
        shape = "rod"
    elif ratios[0] < 0.15:
        shape = "slab"
    else:
        shape = "irregular"
    return rg, asphericity, acylindricity, float(kappa2), shape


# --------------------------------------------------------------------------- #
def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dir", required=True)
    p.add_argument("--traj", required=True)
    p.add_argument("--psf")
    p.add_argument("--last-ns", type=float)
    p.add_argument("--frame-ps", type=float, help="frame spacing for VMD exports")
    p.add_argument("--stride", type=int, default=1)
    p.add_argument("--max-frames", type=int, default=0)
    p.add_argument("--outdir", default="analysis")
    p.add_argument("--cutoff", type=float, default=4.0, help="contact / centroid cutoff in A")
    p.add_argument("--mode", choices=("contact", "centroid"), default="contact")
    p.add_argument("--shell", type=float, default=0.5, help="radial density shell width in A")
    p.add_argument("--rmax", type=float, default=30.0, help="radial density max distance in A")
    p.add_argument("--no-density", action="store_true")
    return p.parse_args(argv)


def main(argv=None):
    a = parse_args(argv)
    d = os.path.abspath(a.dir)
    traj = a.traj if os.path.isabs(a.traj) else os.path.join(d, a.traj)
    if not os.path.exists(traj):
        sys.exit(f"no such trajectory: {traj}")

    case = scaffold.load_case(d) if os.path.exists(os.path.join(d, scaffold.CASE_FILE)) else None
    ts_fs = (case or {}).get("timestep_fs", 2.0)
    psf = a.psf or os.path.join(d, (case or {}).get("psf", "ionized.psf"))

    natoms, atoms = namdlib.read_psf(psf)
    molecules = solute_molecules(atoms)
    if not molecules:
        sys.exit("no solute molecules found (everything looks like water/ions)")
    kinds = Counter(rn for (rn, _rid), _ats in molecules)

    dcd = namdlib.Dcd(traj)
    idx, per, how = index_map(dcd, natoms, molecules)
    nmol = len(molecules)

    dt_ps = ts_fs / 1000.0
    if a.frame_ps is not None:
        ps_per_frame = float(a.frame_ps)
    elif dcd.nsavc and dcd.nsavc > 1:
        ps_per_frame = dcd.nsavc * dt_ps
    else:
        ps_per_frame = 0.0
    ns_per_frame = ps_per_frame / 1000.0

    start = 0
    last = a.last_ns if a.last_ns is not None else (case or {}).get("analysis", {}).get("default_last_ns")
    if last and ns_per_frame > 0:
        start = max(0, dcd.frames - int(round(last / ns_per_frame)))
    frames = list(range(start, dcd.frames, max(1, a.stride)))
    if a.max_frames:
        frames = frames[: a.max_frames]
    if not frames:
        sys.exit("no frames selected")

    print(f"trajectory : {os.path.basename(traj)}")
    print(f"atoms      : {dcd.natoms:,} ({how})")
    print(f"solute     : {nmol} molecules, types {dict(kinds)}")
    span = (f"{start * ns_per_frame:.1f} -> {dcd.frames * ns_per_frame:.1f} ns, "
            f"{ps_per_frame:.1f} ps apart") if ps_per_frame > 0 else "step info unavailable"
    print(f"frames     : {len(frames):,} of {dcd.frames:,}  ({span})")

    outdir = a.outdir if os.path.isabs(a.outdir) else os.path.join(d, a.outdir)
    try:
        os.makedirs(outdir, exist_ok=True)
        probe = os.path.join(outdir, ".write_test")
        open(probe, "w").close()
        os.remove(probe)
    except OSError as exc:
        sys.exit(f"cannot write to {outdir}: {exc}")

    multi = len(kinds) > 1
    size_counts = Counter()
    series = []
    shape_rows = []
    comp_rows = []
    shell_edges = np.arange(0.0, a.rmax + a.shell * 0.5, a.shell)
    dens = np.zeros(len(shell_edges) - 1)
    dens_frames = 0
    charge_state = Counter()

    offsets = np.concatenate([[0], np.cumsum(per)])

    for k, i in enumerate(frames):
        coords, cell = dcd.frame(i, idx)
        if cell is None:
            sys.exit("trajectory has no unit cell; pass a trajectory that carries it")
        x = np.asarray(coords, dtype=np.float64)
        mol_xyz = [x[offsets[m]:offsets[m + 1]] for m in range(nmol)]
        mol_xyz = [namdlib.unwrap_molecule(c, cell) for c in mol_xyz]

        labels = cluster_labels(mol_xyz, cell, a.cutoff, a.mode)
        groups = group_by_label(labels)
        sizes = [len(g) for g in groups]
        for s in sizes:
            size_counts[s] += 1

        t_val = i * ns_per_frame if ps_per_frame > 0 else float(i)
        series.append((t_val, len(groups), max(sizes), float(np.mean(sizes))))

        # largest cluster: shape + composition + radial density
        big = groups[0]
        pts = gather_cluster(mol_xyz, big, cell)
        rg, asp, acyl, kappa2, shape = gyration(pts)
        shape_rows.append((t_val, len(big), rg, asp, acyl, kappa2, shape))

        if multi:
            counts = Counter(molecules[m][0][0] for m in big)
            comp_rows.append((t_val, len(big), dict(counts)))

        if not a.no_density and len(big) >= 2:
            centre = pts.mean(axis=0)
            delta = pts - centre
            delta = namdlib.minimum_image(delta, cell)
            dist = np.linalg.norm(delta, axis=1)
            dens += np.histogram(dist, bins=shell_edges)[0]
            dens_frames += 1

        for m in big:
            charge_state[molecules[m][0][0]] += 1

        if (k + 1) % 200 == 0:
            print(f"  ... {k + 1}/{len(frames)} frames", flush=True)

    dcd.close()

    # ---- outputs ---------------------------------------------------------
    total_frames = len(frames)
    with open(os.path.join(outdir, "clusters.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["cluster_size", "n_clusters", "fraction_of_frames"])
        for s in sorted(size_counts):
            w.writerow([s, size_counts[s], f"{size_counts[s] / total_frames:.5f}"])

    with open(os.path.join(outdir, "clusters_timeseries.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["time_ns" if ps_per_frame > 0 else "frame",
                    "n_clusters", "largest_cluster", "mean_cluster_size"])
        for row in series:
            w.writerow([f"{row[0]:.3f}", row[1], row[2], f"{row[3]:.2f}"])

    with open(os.path.join(outdir, "shape.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["time_ns" if ps_per_frame > 0 else "frame", "largest_cluster",
                    "Rg_A", "asphericity", "acylindricity", "rel_anisotropy", "shape"])
        for row in shape_rows:
            w.writerow([f"{row[0]:.3f}", row[1], f"{row[2]:.2f}", f"{row[3]:.4f}",
                        f"{row[4]:.4f}", f"{row[5]:.4f}", row[6]])

    if comp_rows:
        all_types = sorted(kinds)
        with open(os.path.join(outdir, "composition.csv"), "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["time_ns" if ps_per_frame > 0 else "frame", "largest_cluster"]
                       + [f"n_{t}" for t in all_types])
            for t_val, size, counts in comp_rows:
                w.writerow([f"{t_val:.3f}", size] + [counts.get(t, 0) for t in all_types])

    if dens_frames:
        mid = 0.5 * (shell_edges[:-1] + shell_edges[1:])
        shell_vol = (4.0 / 3.0) * np.pi * (shell_edges[1:] ** 3 - shell_edges[:-1] ** 3)
        density = dens / (dens_frames * shell_vol)
        with open(os.path.join(outdir, "density_profile.csv"), "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["r_A", "solute_density_per_A3"])
            for r_, rho in zip(mid, density):
                w.writerow([f"{r_:.2f}", f"{rho:.6f}"])

    # ---- summary ---------------------------------------------------------
    sizes_flat = [r[2] for r in series]        # r[2] is the largest cluster size
    rgs = [r[2] for r in shape_rows]
    shapes = Counter(r[6] for r in shape_rows)
    q = max(1, len(series) // 4)
    first_mean = np.mean([r[3] for r in series[:q]])
    last_mean = np.mean([r[3] for r in series[-q:]])
    trend = ("growing" if last_mean > first_mean * 1.15
             else "shrinking" if last_mean < first_mean * 0.85 else "stable")

    lines = [
        f"trajectory : {os.path.basename(traj)}",
        f"frames     : {total_frames:,} (stride {a.stride})",
        f"solute     : {nmol} molecules of {dict(kinds)}",
        f"criterion  : {a.mode}, cutoff {a.cutoff} A",
        "",
        f"largest cluster : mean {np.mean(sizes_flat):.1f}, max {max(sizes_flat)}, "
        f"final {sizes_flat[-1]} of {nmol} molecules",
        f"mean cluster    : {first_mean:.1f} -> {last_mean:.1f} molecules ({trend})",
        f"radius of gyration of the largest cluster: mean {np.mean(rgs):.2f} A, "
        f"last {rgs[-1]:.2f} A",
        f"shape of the largest cluster: " + ", ".join(f"{k} {v}" for k, v in shapes.most_common()),
    ]
    if dens_frames:
        peak = float(mid[int(np.argmax(density))])
        lines.append(f"radial density peaks at {peak:.1f} A from the cluster centre "
                     f"(core-shell structure if this is not ~0)")
    if multi:
        lines.append("composition: see composition.csv (molecule types in the largest cluster)")

    txt = "\n".join(lines)
    with open(os.path.join(outdir, "summary.txt"), "w", encoding="utf-8") as fh:
        fh.write(txt + "\n")
    print("\n" + txt)
    print(f"\nwrote {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())