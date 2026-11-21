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

"""Aggregation analysis for a solute-only or full-system NAMD trajectory.

Three measurements, matching the J-aggregation paper's figures:
  * radial distribution functions between named atoms, intermolecular pairs only
  * interplanar angle and distance between two aromatic plane groups
  * aggregate size distribution from a contact criterion

Outputs CSV files (Origin/Grace friendly) plus a text summary.

    python analyze.py --dir <dir> --traj 04_prod_02.dcd --last-ns 50
    python analyze.py --dir <dir> --traj nowater.dcd --planes
    python analyze.py --dir <dir> --traj 04_prod.dcd --stride 10 --no-rdf
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import namdlib  # noqa: E402
import scaffold  # noqa: E402


# --------------------------------------------------------------------------- #
# molecule bookkeeping
# --------------------------------------------------------------------------- #
def load_molecules(psf_path, resname):
    """Return (natoms, {resid: [atom dict...]} ordered by resid)."""
    natoms, atoms = namdlib.read_psf(psf_path)
    mols = defaultdict(list)
    for a in atoms:
        if a["resname"] == resname:
            mols[a["resid"]].append(a)
    mols = {k: mols[k] for k in sorted(mols)}
    if not mols:
        raise SystemExit(f"no residue named {resname} in {psf_path}")
    return natoms, mols


def resolve_indices(dcd, psf_natoms, mols):
    """Map PSF ligand atoms onto trajectory atom slots."""
    n_per = {len(v) for v in mols.values()}
    if len(n_per) != 1:
        raise SystemExit(f"ligand residues have inconsistent atom counts: {sorted(n_per)}")
    per = n_per.pop()
    nmol = len(mols)

    if dcd.natoms == psf_natoms:
        order = []
        for resid, atoms in mols.items():
            order.extend(a["serial"] - 1 for a in atoms)
        return np.array(order, dtype=np.int64), nmol, per, "full system (indices from PSF)"

    if dcd.natoms == nmol * per:
        return np.arange(nmol * per, dtype=np.int64), nmol, per, "solute-only trajectory"

    raise SystemExit(
        f"trajectory has {dcd.natoms:,} atoms; PSF has {psf_natoms:,} and "
        f"{nmol} x {per} = {nmol * per:,} ligand atoms. Neither matches. "
        "Pass the matching PSF, or pre-select with VMD and use that .pdb ordering."
    )


# --------------------------------------------------------------------------- #
# measurements
# --------------------------------------------------------------------------- #
def frame_bytes_for(natoms):  # pragma: no cover - kept for API compatibility
    return namdlib.frame_bytes_for(natoms)


def safe_name(text):
    """Make a site/group name usable as a Windows file name ('ring:a' -> 'ring-a')."""
    return re.sub(r"[^A-Za-z0-9_.-]+", "-", text).strip("-")


def rdf_grid(rmax, dr):
    edges = np.arange(0.0, rmax + dr * 0.5, dr)
    return edges, 0.5 * (edges[:-1] + edges[1:])


def find_peaks(r, g, rmin=2.5, n=3):
    """Local maxima of g(r) above rmin, strongest first, returned sorted by r."""
    out = []
    for i in range(1, len(g) - 1):
        if r[i] < rmin:
            continue
        if g[i] > g[i - 1] and g[i] >= g[i + 1]:
            out.append((float(r[i]), float(g[i])))
    out.sort(key=lambda t: -t[1])
    return sorted(out[:n])


def accumulate_rdf(pos_a, pos_b, cell, hist, edges):
    """Add intermolecular pair distances between two (N,3) arrays.

    Only i<j is counted; the ordered-pair factor of 2 is applied at
    normalisation time. Pairs belonging to the same molecule are always
    excluded, which is what makes the RDF purely intermolecular.
    """
    n = len(pos_a)
    delta = pos_a[:, None, :] - pos_b[None, :, :]
    delta = namdlib.minimum_image(delta, cell)
    dist = np.linalg.norm(delta, axis=2)
    iu = np.triu_indices(n, k=1)
    h, _ = np.histogram(dist[iu], bins=edges)
    hist += h


def plane_frame(coords_by_mol, group_atom_idx, cell):
    """Per-molecule centroid, normal, planarity rms for one atom group."""
    cents, normals, rms = [], [], []
    for x in coords_by_mol:
        pts = x[group_atom_idx]
        c, n, r = namdlib.plane(pts)
        cents.append(c)
        normals.append(n)
        rms.append(r)
    return np.array(cents), np.array(normals), np.array(rms)


def cluster_sizes(centroids, coords_by_mol, cell, cutoff, mode="contact", pad=12.0):
    """Union-find over molecules; returns list of component sizes."""
    n = len(centroids)
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

    d = centroids[:, None, :] - centroids[None, :, :]
    d = namdlib.minimum_image(d, cell)
    cen_dist = np.linalg.norm(d, axis=2)
    rows, cols = np.where(np.triu(cen_dist < cutoff + pad, k=1))
    for i, j in zip(rows, cols):
        if mode == "centroid":
            if cen_dist[i, j] <= cutoff:
                union(i, j)
            continue
        delta = coords_by_mol[i][:, None, :] - coords_by_mol[j][None, :, :]
        delta = namdlib.minimum_image(delta, cell)
        if np.sqrt((delta ** 2).sum(axis=2)).min() <= cutoff:
            union(i, j)

    sizes = defaultdict(int)
    for i in range(n):
        sizes[find(i)] += 1
    return sorted(sizes.values(), reverse=True)


# --------------------------------------------------------------------------- #
def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dir", required=True)
    p.add_argument("--traj", required=True, help=".dcd file name inside --dir")
    p.add_argument("--psf", help="PSF to map atoms (default: case.json / ionized.psf)")
    p.add_argument("--last-ns", type=float, help="only the last N ns (default: case.json value)")
    p.add_argument("--stride", type=int, default=1)
    p.add_argument("--max-frames", type=int, default=0)
    p.add_argument("--outdir", default="analysis",
                   help="output directory; relative paths are inside --dir, absolute paths win. "
                        "Point this somewhere writable if the case directory is read-only.")
    p.add_argument("--frame-ps", type=float,
                   help="frame spacing for trajectories that carry no step information "
                        "(e.g. a .dcd written by VMD)")
    p.add_argument("--rmax", type=float, default=12.0)
    p.add_argument("--dr", type=float, default=0.05)
    p.add_argument("--rdf", help="pairs like C5:O3,C5:C5 (default from case.json)")
    p.add_argument("--planes", nargs=2, metavar=("GROUP_A", "GROUP_B"),
                   help="plane groups to compare (default from case.json)")
    p.add_argument("--plane-mode", choices=("inter", "intra"), default="inter",
                   help="inter = planes of different molecules (J-aggregate geometry, default); "
                        "intra = the two groups inside one molecule")
    p.add_argument("--plane-cutoff", type=float, default=12.0,
                   help="inter mode: only molecules whose centroids are closer than this (A)")
    p.add_argument("--slip-threshold", type=float, default=32.0,
                   help="slip angle (deg) separating H-type (<) from J-type (>); default 32")
    p.add_argument("--lateral-cutoff", type=float, default=6.0,
                   help="max in-plane offset (A) for a pair to count as a stack rather than "
                        "two neighbours lying side by side in the same layer")
    p.add_argument("--stack-cutoff", type=float, default=8.0,
                   help="max centre-to-centre distance (A) for a pair to count as adjacent "
                        "stacked molecules")
    p.add_argument("--cluster-cutoff", type=float, help="contact cutoff in A")
    p.add_argument("--cluster-mode", choices=("contact", "centroid"), default="contact")
    p.add_argument("--box", help="box lengths a,b,c for trajectories without a unit-cell record")
    p.add_argument("--no-rdf", action="store_true")
    p.add_argument("--no-planes", action="store_true")
    p.add_argument("--no-clusters", action="store_true")
    return p.parse_args(argv)


def main(argv=None):
    a = parse_args(argv)
    d = os.path.abspath(a.dir)
    traj = a.traj if os.path.isabs(a.traj) else os.path.join(d, a.traj)
    if not os.path.exists(traj):
        sys.exit(f"no such trajectory: {traj}")

    case = scaffold.load_case(d) if os.path.exists(os.path.join(d, scaffold.CASE_FILE)) else None
    ana = (case or scaffold.DEFAULT_CASE)["analysis"]
    resname = (case or {}).get("ligand_resname", "LIG")
    ts_fs = (case or {}).get("timestep_fs", 2.0)
    psf = a.psf or os.path.join(d, (case or {}).get("psf", "ionized.psf"))

    psf_natoms, mols = load_molecules(psf, resname)
    dcd = namdlib.Dcd(traj)
    indices, nmol, per, how = resolve_indices(dcd, psf_natoms, mols)
    dt_ps = ts_fs / 1000.0
    # A DCD written by VMD carries nsavc=1 (one MD step), which is not a frame spacing.
    # Only trust it above one step, or take --frame-ps as the explicit answer.
    if a.frame_ps is not None:
        ps_per_frame = float(a.frame_ps)
    elif dcd.nsavc and dcd.nsavc > 1:
        ps_per_frame = dcd.nsavc * dt_ps
    else:
        ps_per_frame = 0.0
    ns_per_frame = ps_per_frame / 1000.0

    last_ns = a.last_ns if a.last_ns is not None else ana.get("default_last_ns")
    start = 0
    if last_ns:
        if ps_per_frame <= 0:
            print("[warn] this trajectory has no step information (nsavc=0), so --last-ns "
                  "cannot be applied; analysing every selected frame. Pass --frame-ps to declare "
                  "the spacing.")
            last_ns = None
        else:
            n_back = int(round(last_ns / ns_per_frame))
            start = max(0, dcd.frames - n_back)
    frames = list(range(start, dcd.frames, max(1, a.stride)))
    if a.max_frames:
        frames = frames[: a.max_frames]

    if ps_per_frame > 0:
        t0 = start * ns_per_frame
        t1 = dcd.frames * ns_per_frame
        span = f"{t0:.1f} -> {t1:.1f} ns, {ps_per_frame:.1f} ps apart"
    else:
        span = "step information unavailable (VMD export)"
    print(f"trajectory : {os.path.basename(traj)}")
    print(f"atoms      : {dcd.natoms:,} ({how})")
    print(f"molecules  : {nmol} x {per} atoms ({resname})")
    print(f"frames     : {len(frames):,} of {dcd.frames:,}  ({span})")
    if not frames:
        sys.exit("no frames selected")

    fallback_cell = None
    if not dcd.has_cell:
        if a.box:
            lens = [float(x) for x in a.box.replace(",", " ").split()]
            if len(lens) != 3:
                sys.exit("--box needs three numbers: a,b,c")
            fallback_cell = np.diag(lens)
        else:
            pdb_path = os.path.join(d, (case or {}).get("pdb", "ionized.pdb"))
            cryst = namdlib.read_pdb_box(pdb_path) if os.path.exists(pdb_path) else None
            if cryst:
                fallback_cell = np.diag([cryst[0], cryst[1], cryst[2]])
                print(f"box        : no unit cell in the trajectory, using CRYST1 from "
                      f"{os.path.basename(pdb_path)} ({cryst[0]:g} {cryst[1]:g} {cryst[2]:g})")
            else:
                sys.exit("this trajectory has no unit-cell record and no CRYST1 was found; "
                         "pass --box a,b,c")

    outdir = a.outdir if os.path.isabs(a.outdir) else os.path.join(d, a.outdir)
    try:
        os.makedirs(outdir, exist_ok=True)
        probe = os.path.join(outdir, ".write_test")
        with open(probe, "w", encoding="utf-8") as fh:
            fh.write("x")
        os.remove(probe)
    except OSError as exc:
        sys.exit(f"cannot write to {outdir}: {exc}\n"
                 f"Pick a writable location, e.g. --outdir {os.path.join(os.path.dirname(d), '_analysis', os.path.basename(d))}")

    # molecule -> [atom names]
    mol_atom_names = {k: [x["name"] for x in v] for k, v in mols.items()}
    first_mol = next(iter(mol_atom_names))
    name_to_local = {n: i for i, n in enumerate(mol_atom_names[first_mol])}

    # --- prepare requests ---------------------------------------------------
    raw_groups = ana.get("plane_groups", {})

    def split_pair(token):
        """'C5|O3', legacy 'C5:O3', or 'ring:ring1' (centre of a group)."""
        token = token.strip()
        if "|" in token:
            left, right = token.split("|", 1)
            return left.strip(), right.strip()
        if token.startswith("ring:") and token.count(":") == 1:
            return token, token
        parts = token.split(":")
        if len(parts) == 2:
            return parts[0].strip(), parts[1].strip()
        raise ValueError(f"cannot read RDF pair {token!r}; use 'A|B', 'A:B' or 'ring:NAME'")

    def site_ok(site):
        if site.startswith("ring:"):
            return site[5:] in raw_groups
        return site in name_to_local

    pairs, rdf_sites = [], set()
    if not a.no_rdf:
        spec = ana.get("rdf_pairs", [])
        if a.rdf:
            spec = [split_pair(t) for t in a.rdf.split(",") if t.strip()]
        for x, y in spec:
            if not (site_ok(x) and site_ok(y)):
                missing = [s for s in (x, y) if not site_ok(s)]
                print(f"[skip] rdf {x}|{y} -- unknown site(s) {missing}; "
                      "atom names must exist in the solute, ring groups must be in case.json")
                continue
            pairs.append((x, y))
            rdf_sites.update((x, y))

    groups = {}
    plane_pair = None
    if not a.no_planes:
        if a.planes:
            plane_pair = tuple(a.planes)
        else:
            plane_pair = tuple(ana.get("plane_pair", [])) or None
        if plane_pair:
            for g in plane_pair:
                names = raw_groups.get(g)
                if not names:
                    print(f"[skip] plane group '{g}' is not defined in case.json")
                    plane_pair = None
                    break
                missing = [n for n in names if n not in name_to_local]
                if missing:
                    sys.exit(f"plane group '{g}' uses atom names not in the ligand: {missing}")
                groups[g] = np.array([name_to_local[n] for n in names], dtype=np.int64)

    # ring sites used by the RDF but not requested as planes
    for site in sorted(rdf_sites):
        if site.startswith("ring:"):
            g = site[5:]
            if g not in groups:
                groups[g] = np.array([name_to_local[n] for n in raw_groups[g]], dtype=np.int64)

    if not pairs and not plane_pair and not raw_groups:
        print("\n[warn] nothing to analyse: case.json has no plane_groups / rdf_pairs.\n"
              "       Run  inspect_molecule.py --dir <dir> --write  to detect ring systems and\n"
              "       fill them in, or pass --planes A B / --rdf 'A|B' explicitly.\n"
              "       Cluster analysis still runs.")

    cutoff = a.cluster_cutoff if a.cluster_cutoff else ana.get("cluster_cutoff_A", 4.0)

    # --- accumulators -------------------------------------------------------
    edges, centers = rdf_grid(a.rmax, a.dr)
    hist = {p: np.zeros(len(centers)) for p in pairs}
    angle_hist = np.zeros(91)
    slip_hist = np.zeros(91)
    dist_edges = np.arange(0.0, 20.0001, 0.25)
    dist_hist = np.zeros(len(dist_edges) - 1)
    nearest_edges = np.arange(0.0, 10.0001, 0.1)
    nearest_hist = np.zeros(len(nearest_edges) - 1)
    stack_hist = np.zeros(len(dist_edges) - 1)
    stack_slip_sum = 0.0
    stack_n = 0
    cluster_counts = defaultdict(int)
    cluster_series = []
    volume_sum = 0.0
    n_used = 0
    plane_rms = []
    n_planes = 0

    # --- loop ---------------------------------------------------------------
    for k, i in enumerate(frames):
        coords, cell = dcd.frame(i, indices)
        if cell is None:
            cell = fallback_cell
        x = np.asarray(coords, dtype=np.float64)
        vol = abs(np.linalg.det(cell))
        volume_sum += vol
        coords_by_mol = [x[m * per:(m + 1) * per] for m in range(nmol)]

        # unwrap each molecule for centroids / planes
        unc = [namdlib.unwrap_molecule(c, cell) for c in coords_by_mol]
        centroids = np.array([u.mean(axis=0) for u in unc])

        if pairs:
            site_pos = {}
            for site in rdf_sites:
                if site.startswith("ring:"):
                    gidx = groups[site[5:]]
                    site_pos[site] = np.array([u[gidx].mean(axis=0) for u in unc])
                else:
                    idxs = np.array([name_to_local[site] + m * per for m in range(nmol)])
                    site_pos[site] = x[idxs]
            for key in pairs:
                accumulate_rdf(site_pos[key[0]], site_pos[key[1]], cell, hist[key], edges)

        if plane_pair:
            cents_a, normals_a, rms_a = plane_frame(unc, groups[plane_pair[0]], cell)
            cents_b, normals_b, rms_b = plane_frame(unc, groups[plane_pair[1]], cell)
            plane_rms.append(float(np.concatenate([rms_a, rms_b]).mean()))
            if a.plane_mode == "intra":
                ii = np.arange(nmol)
                jj = ii.copy()
            else:
                dcen = cents_a[:, None, :] - cents_b[None, :, :]
                dcen = namdlib.minimum_image(dcen, cell)
                cen_dist = np.linalg.norm(dcen, axis=2)
                mask = np.triu(cen_dist < a.plane_cutoff, k=1)
                ii, jj = np.where(mask)
            if len(ii):
                cosang = np.abs(np.clip(
                    np.einsum("ij,ij->i", normals_a[ii], normals_b[jj]), -1.0, 1.0))
                ang = np.degrees(np.arccos(cosang))
                angle_hist += np.bincount(np.clip(ang.astype(int), 0, 90), minlength=91)[:91]
                dvec = cents_b[jj] - cents_a[ii]
                dvec = namdlib.minimum_image(dvec, cell)
                dist = np.abs(np.einsum("ij,ij->i", dvec, normals_a[ii]))
                dist_hist += np.histogram(dist, bins=dist_edges)[0]
                # Slip angle: angle between the centroid-centroid vector and the plane
                # normal. Near 0 means face-to-face (H-type); larger means the molecules
                # are slipped along the plane (J-type). This, not the normal-normal
                # angle, is what separates the two packing modes.
                cdist = np.linalg.norm(dvec, axis=1)
                proj = np.abs(np.einsum("ij,ij->i", dvec, normals_a[ii]))
                slip = np.degrees(np.arccos(np.clip(proj / np.maximum(cdist, 1e-9), 0.0, 1.0)))
                slip_hist += np.bincount(np.clip(slip.astype(int), 0, 90), minlength=91)[:91]
                # A small perpendicular offset alone does NOT mean two molecules are
                # stacked: neighbours lying side by side in the same layer also have a
                # near-zero perpendicular component, and molecules far apart along the
                # normal have a large one. A stack needs BOTH a bounded centre distance
                # and a bounded lateral offset.
                lateral = np.sqrt(np.maximum(cdist ** 2 - dist ** 2, 0.0))
                stacked = (lateral <= a.lateral_cutoff) & (cdist <= a.stack_cutoff)
                if stacked.any():
                    stack_hist += np.histogram(dist[stacked], bins=dist_edges)[0]
                    stack_slip_sum += float(slip[stacked].sum())
                    stack_n += int(stacked.sum())
                # Spacing to the nearest neighbouring molecule (centre-to-centre),
                # which is what "distance between adjacent molecules" means.
                per_mol = np.full(nmol, np.inf)
                np.minimum.at(per_mol, ii, cdist)
                np.minimum.at(per_mol, jj, cdist)
                ok = np.isfinite(per_mol)
                nearest_hist += np.histogram(per_mol[ok], bins=nearest_edges)[0]
                n_planes += len(ii)

        if not a.no_clusters:
            sizes = cluster_sizes(centroids, unc, cell, cutoff, a.cluster_mode)
            for size in sizes:
                cluster_counts[size] += 1
            t_val = (i * ns_per_frame) if ps_per_frame > 0 else float(i)
            cluster_series.append((t_val, len(sizes),
                                   max(sizes) if sizes else 0,
                                   (sum(sizes) / len(sizes)) if sizes else 0.0))

        n_used += 1
        if n_used % 200 == 0:
            print(f"  ... {n_used}/{len(frames)} frames", flush=True)

    dcd.close()
    avg_vol = volume_sum / n_used
    summary = []

    # --- RDF output ---------------------------------------------------------
    for (xname, yname), h in hist.items():
        h = h.astype(float)
        pref = 2.0                      # i<j counted once, ordered pairs count twice
        rho = (nmol - 1) / avg_vol      # every molecule contributes one of each atom
        shell = 4.0 * np.pi * centers ** 2 * a.dr
        g = np.where(shell > 0, pref * h / (n_used * nmol * rho * np.where(shell > 0, shell, 1)), 0.0)
        fname = os.path.join(outdir, f"rdf_{safe_name(xname)}_{safe_name(yname)}.csv")
        with open(fname, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["r_A", "g(r)", "count", "frames"])
            for r_, g_, c_ in zip(centers, g, h):
                w.writerow([f"{r_:.3f}", f"{g_:.4f}", int(c_), n_used])
        peaks = find_peaks(centers, g, rmin=2.5)
        desc = ", ".join(f"r={r_:.2f} g={g_:.2f}" for r_, g_ in peaks) or "no peak above 2.5 A"
        summary.append(f"RDF {xname}-{yname}: {desc} ({n_used} frames) -> {os.path.basename(fname)}")

    # --- planes output ------------------------------------------------------
    if plane_pair:
        ga, gb = plane_pair
        ang_mid = np.arange(0.5, 91.0, 1.0)
        f1 = os.path.join(outdir, f"planes_{safe_name(ga)}_{safe_name(gb)}_angle.csv")
        with open(f1, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["angle_deg", "count", "fraction"])
            for m_, c_ in zip(ang_mid, angle_hist):
                w.writerow([f"{m_:.1f}", int(c_), f"{c_ / max(1, angle_hist.sum()):.5f}"])
        dmid = 0.5 * (dist_edges[:-1] + dist_edges[1:])
        f2 = os.path.join(outdir, f"planes_{safe_name(ga)}_{safe_name(gb)}_distance.csv")
        with open(f2, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["d_A", "count", "fraction"])
            for m_, c_ in zip(dmid, dist_hist):
                w.writerow([f"{m_:.3f}", int(c_), f"{c_ / max(1, dist_hist.sum()):.5f}"])

        mean_ang = float((ang_mid * angle_hist).sum() / max(1, angle_hist.sum()))
        mean_d = float((dmid * dist_hist).sum() / max(1, dist_hist.sum()))
        frac_d = float(dist_hist[: int(4.0 / 0.25)].sum() / max(1, dist_hist.sum()))
        slip_mid = np.arange(0.5, 91.0, 1.0)
        mean_slip = float((slip_mid * slip_hist).sum() / max(1, slip_hist.sum()))
        thr = int(round(a.slip_threshold))
        frac_jtype = float(slip_hist[thr:].sum() / max(1, slip_hist.sum()))
        f3 = os.path.join(outdir, f"planes_{safe_name(ga)}_{safe_name(gb)}_slip.csv")
        with open(f3, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["slip_deg", "count", "fraction"])
            for m_, c_ in zip(slip_mid, slip_hist):
                w.writerow([f"{m_:.1f}", int(c_), f"{c_ / max(1, slip_hist.sum()):.5f}"])
        summary.append(
            f"planes {ga} vs {gb} ({a.plane_mode}, {n_planes:,} pairs): "
            f"mean angle {mean_ang:.1f} deg, mean |d| {mean_d:.2f} A "
            f"({frac_d * 100:.1f}% within 4 A), planarity rms {np.mean(plane_rms):.3f} A")
        summary.append(
            f"slip angle: mean {mean_slip:.1f} deg, J-type (> {thr} deg) "
            f"{frac_jtype * 100:.1f}%, H-type {100 - frac_jtype * 100:.1f}% "
            f"-> {os.path.basename(f3)}")
        nmid = 0.5 * (nearest_edges[:-1] + nearest_edges[1:])
        ntot = max(1, nearest_hist.sum())
        cdf = np.cumsum(nearest_hist) / ntot
        med_r = float(nmid[np.searchsorted(cdf, 0.5)]) if ntot > 1 else float("nan")
        f5 = os.path.join(outdir, f"planes_{safe_name(ga)}_{safe_name(gb)}_nearest.csv")
        with open(f5, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["nearest_centre_distance_A", "count", "fraction"])
            for m_, c_ in zip(nmid, nearest_hist):
                w.writerow([f"{m_:.2f}", int(c_), f"{c_ / ntot:.5f}"])
        summary.append(
            f"nearest neighbour: median {med_r:.2f} A, mean "
            f"{float((nmid * nearest_hist).sum() / ntot):.2f} A between molecular centres "
            f"-> {os.path.basename(f5)}")
        if med_r < 2.0:
            summary.append(
                "SANITY FAIL: a nearest-neighbour distance below 2 A is impossible; the "
                "coordinates are almost certainly not in a common frame (check "
                "namdlib.unwrap_molecule) or the atom selection is wrong")

        dmid2 = 0.5 * (dist_edges[:-1] + dist_edges[1:])
        if stack_n:
            sd_mean = float((dmid2 * stack_hist).sum() / stack_n)
            cdf2 = np.cumsum(stack_hist) / stack_n
            sd_med = float(dmid2[np.searchsorted(cdf2, 0.5)])
            sd_frac4 = float(stack_hist[dmid2 <= 4.0].sum() / stack_n) * 100
            sd_slip = stack_slip_sum / stack_n
            summary.append(
                f"stacked pairs (lateral <= {a.lateral_cutoff:g} A, {stack_n:,}): "
                f"interplanar d median {sd_med:.2f} A, mean {sd_mean:.2f} A, "
                f"{sd_frac4:.1f}% within 4 A, mean slip {sd_slip:.1f} deg")
            f6 = os.path.join(outdir, f"planes_{safe_name(ga)}_{safe_name(gb)}_stacked_d.csv")
            with open(f6, "w", newline="", encoding="utf-8") as fh:
                w = csv.writer(fh)
                w.writerow(["interplanar_d_A", "count", "fraction"])
                for m_, c_ in zip(dmid2, stack_hist):
                    w.writerow([f"{m_:.3f}", int(c_), f"{c_ / stack_n:.5f}"])
        else:
            summary.append("stacked pairs: none within the lateral cutoff")

    # --- clusters output ----------------------------------------------------
    if not a.no_clusters:
        total = sum(cluster_counts.values())
        f3 = os.path.join(outdir, "clusters.csv")
        with open(f3, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["cluster_size", "n_clusters", "fraction_of_frames"])
            for size in sorted(cluster_counts):
                w.writerow([size, cluster_counts[size], f"{cluster_counts[size] / total:.5f}"])
        big = sum(c for s, c in cluster_counts.items() if s >= 2)
        summary.append(f"clusters (cutoff {cutoff} A, {a.cluster_mode}): "
                       f"{big / max(1, total) * 100:.1f}% of frames have clusters, "
                       f"largest ever {max(cluster_counts) if cluster_counts else 0} molecules "
                       f"-> clusters.csv")
        f4 = os.path.join(outdir, "clusters_timeseries.csv")
        with open(f4, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["time_ns" if ps_per_frame > 0 else "frame",
                        "n_clusters", "largest_cluster", "mean_cluster_size"])
            for t, nc, lg, mean in cluster_series:
                w.writerow([f"{t:.3f}", nc, lg, f"{mean:.2f}"])
        if len(cluster_series) >= 4:
            q = max(1, len(cluster_series) // 4)
            first = np.mean([c[3] for c in cluster_series[:q]])
            last = np.mean([c[3] for c in cluster_series[-q:]])
            trend = "growing" if last > first * 1.15 else (
                "shrinking" if last < first * 0.85 else "stable")
            summary.append(f"assembly trend: mean cluster {first:.1f} -> {last:.1f} molecules "
                           f"({trend}) -> {os.path.basename(f4)}")

    txt = "\n".join([
        f"trajectory : {os.path.basename(traj)}",
        f"frames     : {n_used:,} (stride {a.stride})",
        f"box volume : {avg_vol:.0f} A^3 (mean)",
        "",
        *summary,
    ])
    with open(os.path.join(outdir, "summary.txt"), "w", encoding="utf-8") as fh:
        fh.write(txt + "\n")
    print("\n" + txt)
    print(f"\nwrote {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
