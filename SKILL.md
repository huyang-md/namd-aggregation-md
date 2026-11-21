---
name: namd-aggregation-md
description: 'Use for NAMD aggregation MD of small molecules (probe dyes and similar solutes packed with Packmol into water) inside a J-aggregate / photothermal study: scaffold consistent 01_min/02_nvt/03_equil_npt/04_prod .conf files from one case.json, generate Packmol and VMD solvate inputs at a fixed concentration, pre-flight a run before submitting it, build a correct continuation .conf after an interruption, audit .dcd trajectories against their logs, check that sibling systems are actually comparable, decide whether the study can support each claim yet, and compute RDF / J-versus-H slip geometry / interplanar distance / aggregate sizes. Also use when asked how far a run got, why a run died, whether two systems can be compared, or whether an "X aggregates more than Y" statement is supported.'
---

# NAMD aggregation MD

Scripts live in `scripts/`, next to this file. Run them with any Python 3.10+ that has
numpy (in DSH, take `python` from `load_workspace_dependencies`). Nothing else is
required: the PSF/PDB/DCD/XSC readers are plain Python.

```powershell
$py = '<python from load_workspace_dependencies>'
$S  = '<this skill directory>\scripts'
```

**`reference/protocol.md` holds the science**: the claim chain of a
J-aggregate/photothermal study, which deliverable MD is responsible for, the numeric
thresholds behind each claim, and the five gates a substituent comparison must pass.
Read it before interpreting any number this skill prints.

**`README.md` is the human quick start** — installation, a five-minute walkthrough, what
to change for a new molecule, and the three self-checks. Point users there when they ask
how to use the skill rather than how a number was produced.

**Two analysis layers, and the split matters.** `aggregate.py` (command 5a) is the generic
one: cluster sizes, shape, composition, radial density, valid for micelles, vesicles,
fibrils, co-assemblies and inorganic clusters. `analyze.py` (command 5b) is the
specialised one: interplanar distances, slip angles, J-versus-H classification — only
meaningful when the solute has planar aromatic groups. When asked about a non-chromophore
system, use 5a and do not force the stacking machinery onto it.

**`docs/tutorial.md` is the end-to-end reproducible protocol** — 13 steps from deciding
the charge state through to the claim gate, each with an acceptance criterion and a list
of the errors that actually happen. Load it when someone is starting a system from
scratch, writing methods for a paper, or handing the work to another person. It is the
only place the non-script steps (CHARMM-GUI ligand reader, Packmol, psfgen, VMD solvate)
are written down with their verification.

**`docs/checklist.md` is the tickable version** — hand it to whoever is running the work,
and point at it whenever the question is "what do I still have to do myself". The `👤`
items are the four places a human is genuinely required.

**Adopting an existing project** (hand-written `.conf` files, already-run trajectories)
needs no rebuild: `scaffold.py --dir <dir> --case-only` writes `case.json` alone and
leaves every `.conf` untouched, after which `preflight.py`, `compare_systems.py` and
`study_report.py` all work. Drift checking only applies to files this skill generated, so
a hand-written `.conf` is not reported as drifted.

## The assembly route this skill assumes

One molecule at a time, not a protein with a ligand in a pocket:

1. single molecule → CHARMM-GUI Ligand Reader (or ParamChem) → CGenFF `.str`/`.prm`;
2. **Packmol** (`wa.inp`: `structure <mol>.pdb / number 50 / inside box 0 0 0 120 120 120`)
   → `<case>_com.pdb`;
3. `psfgen` over the packed pdb → `<case>_proc.pdb`/`.psf` (50 × atoms-per-molecule);
4. **VMD `solvate`** → `<case>_solv.pdb`/`.psf`;
5. ions → `ionized.pdb`/`.psf`;
6. NAMD `01_min` → `02_nvt` → `03_equil_npt` → `04_prod`.

Steps 1–5 are where the systems silently stop being comparable, so this skill controls
the box in step 2/4 and audits the result in step 6.

## Why this exists

A multi-system project (five substituent variants, 50 molecules each, CHARMM36/CGenFF,
NAMD) ended up with eleven hand-copied `.conf` files whose force-field block was identical
but whose knobs had silently drifted apart. The failure modes below are all taken from
that project:

- **concentrations differed by more than 1.5×** because Packmol's output box is
  the *bounding box of the packed coordinates*, not the box you asked for, and VMD then
  pads around whatever it got. Aggregate size tracked concentration as much as chemistry;
- **`1-4scaling` split the project in two**: three systems ran production at 1.0, two at
  **1.07**, and both of those switched to 1.07 only for production, after equilibrating at
  1.0 — a step change in the potential at the handoff;
- a restart conf recorded a step that did not match its own `.xsc` (nothing detects this
  for you unless you read the `.xsc`);
- a `.conf` named after one stage wrote `outputName` belonging to another;
- frame spacing changed from 100 ps to 10 ps mid-project, making systems incomparable;
- a 100 ns run died partway with `FATAL ERROR: No space left on device`;
- `.BAK` files that looked like backups were 10 % copies of the real trajectory.

Every one of these is checked or prevented by a script here. Do not hand-edit a generated
`.conf`: change `case.json` and re-run `scaffold.py`.

## The ten commands

### Shortcut: steps ③–⑧ in one call

```powershell
& $py "$S\pipeline.py" --dir <case> --name <name> --ligand-pdb mol.pdb `
      --molecules 50 --target-conc-mM 50 --prod-ns 100
```

Sizes the box from the target concentration, runs Packmol, `psfgen`, VMD
`solvate` + `autoionize`, writes the four `.conf` files, then runs the pre-flight and
reports. It stops with an instruction when a tool is missing instead of half-finishing, and
`--skip-pack` / `--skip-build` / `--skip-solvate` resume from a later stage.

Two steps stay manual by design: **the charge state of the ligand** (a chemistry fact, not
something a script can decide) and **the CGenFF parameters** (a web service or the local
CGenFF program). Everything after them is automated. See `docs/tutorial.md`.

The numbered commands below are the same operations done one at a time, for when you need
to inspect or override an intermediate step.

### 0. Teach it your molecule (do this first on a new system)

```powershell
& $py "$S\inspect_molecule.py" --dir <dir>            # show what it finds
& $py "$S\inspect_molecule.py" --dir <dir> --write    # fill case.json
```

Nothing about a molecule is hard-coded in this skill. `inspect_molecule.py` reads the
topology, finds rings, fuses rings that share an edge into ring systems, measures their
planarity from the coordinates, and writes the usable ones into `case.json` as
`plane_groups`, plus a `plane_pair` (largest group against itself) and an `rdf_pairs`
entry that uses **ring centroids** instead of atom names:

```json
"analysis": {
  "plane_groups": {"ring1": ["C1","C2",...], "ring3": [...]},
  "plane_pair": ["ring1","ring1"],
  "rdf_pairs": [["ring:ring1","ring:ring1"]]
}
```

Ring systems that are not planar (they contain sp3 carbons) are dropped automatically;
`--max-rms` (default 0.2 A) controls that. Check the printed table: a group marked usable
with a large rms means the atom set is wrong, and every angle computed from it is noise.
Adjust `plane_pair` by hand if the comparison you need is group A versus group B rather
than self-stacking.

`analyze.py` accepts `ring:NAME` anywhere an atom name is accepted, so a new molecule
needs no renaming — only this one command.

### 1. Scaffold a case

```powershell
& $py "$S\scaffold.py" --dir <system dir> --name <case> --prod-ns 100 --dcdfreq-ps 10
& $py "$S\scaffold.py" --dir <dir> --dry-run          # print without writing
& $py "$S\scaffold.py" --dir <dir> --force            # adopt the recipe in an existing dir
```

Writes `case.json`, `01_min.conf`, `02_nvt.conf`, `03_equil_npt.conf`, `04_prod.conf`.
It reads `ionized.psf` for the atom count, auto-detects the solute residue (the largest
residue that is neither water nor an ion; override with `--resname`), and takes the box
from `01_min.conf` if present, otherwise from the PDB `CRYST1` (3 decimals — pass
`--cell ax,ay,az,...,ox,oy,oz` when exact box values matter). It also prints the disk the
production run will need.

`case.json` holds everything tunable: parameter list, timestep, `dcdfreq_ps`,
temperature, stage lengths, solute residue name, and the analysis settings. The analysis
block starts **empty on purpose** — fill it with command 0, never by copying atom names
from another molecule.

#### Packing at a fixed concentration

```powershell
& $py "$S\scaffold.py" --dir <new dir> --name <case> --packmol `
      --molecules 50 --target-conc-mM 50 --ligand-pdb <mol>.pdb --dry-run
```

Writes `wa.inp` and `solvate.tcl` as well. The final box edge is computed from the target
concentration, and the Packmol box is deliberately smaller by the water-shell padding, so
the **final** volume — not the bounding box of the packed coordinates — sets the
concentration. Without this, Packmol returns the extent of whatever it placed (50 small
molecules span ~109 A, 50 large ones ~121 A in the same 120 A request) and the systems
end up 1.5x apart in concentration.

Check the generated `solvate.tcl` before running it: it uses `-minmax` to pin the box, and
the PSF argument must be the post-`psfgen` file name for your case.

### 2. Pre-flight before submitting

```powershell
& $py "$S\preflight.py" --dir <dir>            # exit 2 when something is FAIL
& $py "$S\preflight.py" --dir <dir> --verbose  # also show the passing checks
```

Checks, in order of what actually bites: parameter files resolve; PSF/PDB present and
atom counts agree; box vs `CRYST1`; `1-4scaling`/`cutoff`/`timestep`/`exclude` identical
across every `.conf`; `temperature` not combined with `binVelocities` (NAMD aborts);
`cellBasisVector` not combined with `extendedSystem`; `firsttimestep` against the `.xsc`
the conf actually reads; `outputName` not about to append to an existing trajectory;
free disk against the projected trajectory size; and, when `case.json` exists,
whether any `.conf` has drifted from the recipe.

**Read the FAIL lines before launching.** A WARN on `firsttimestep 0` is informational:
a fresh stage legitimately restarts numbering at 0, but remember that a later
continuation of that stage must use the step from its `.xsc`.

### 3. Continuation after an interruption

```powershell
& $py "$S\restart.py" --dir <dir> --from prod_final --target-ns 100               # preview
& $py "$S\restart.py" --dir <dir> --from prod_final --target-ns 100 --write       # save .conf
& $py "$S\restart.py" --dir <dir> --from prod_final --target-ns 100 --part-ns 20 --write
```

Derives `firsttimestep` from the `.xsc` itself, renames the output, drops
`cellBasisVector`, sets `run` to what is left (or to `--part-ns` for queue time limits),
and prints a disk estimate.

It also compares the three restart files' timestamps. NAMD writes `.xsc` *before*
`.coor`/`.vel`, so a crash can leave a triplet from two different steps; when the spread
exceeds two minutes the script warns, because continuing from that triplet restarts from
the `.xsc` step with older geometry. Prefer a triplet written together, or re-run from a
step where all three agree.

### 4. Trajectory audit

```powershell
& $py "$S\inspect_traj.py" --dir <dir>
& $py "$S\inspect_traj.py" --dir <dir> --all      # include equilibration files too
```

Reads each DCD header (NSET / ISTART / NSAVC / NSTEP), counts frames on disk, matches
them against the log's last `ENERGY:` line, then merges the *production* segments per
system and prints covered ns, gaps and overlaps.

Use it whenever someone asks "is this trajectory complete". The header is authoritative:
NAMD rewrites NSET/NSTEP as it writes, so header and disk count agreeing means every
frame landed, and `partial_bytes` at the end means the last frame was cut mid-write.
Never conclude from the `.log` alone — a buffered log can lag the file, and a clean
`End of program` says nothing about how much time was simulated.

Mixed frame spacing inside one system is reported explicitly: averaging frames across
100 ps and 10 ps segments silently gives the coarse segment one tenth of the weight it
should have. Weight by time, not by frame.

### 5a. Aggregation structure — any self-assembly system

```powershell
& $py "$S\aggregate.py" --dir <dir> --traj 04_prod.dcd --last-ns 50
& $py "$S\aggregate.py" --dir <dir> --traj nowater.dcd --frame-ps 100 --mode contact
```

**This is the part that generalises beyond chromophores.** It finds every residue that is
neither water nor an ion and treats it as a solute molecule, so it works unchanged for
micelles and vesicles, fibrils and gels, two-component co-assemblies and inorganic or
hybrid clusters. No `case.json`, no atom names, no ring detection needed.

It reports, per sampled frame: cluster size distribution, the largest and mean cluster,
the **radius of gyration** of the largest cluster, three **shape descriptors**
(asphericity, acylindricity, relative anisotropy) with a sphere/rod/slab/irregular label,
and for multi-component systems the **composition** of each cluster. It also accumulates a
**radial density profile** of solute around the largest cluster's centre — the cheap way to
tell a uniform aggregate from a core-shell micelle or vesicle.

Outputs `clusters.csv`, `clusters_timeseries.csv`, `shape.csv`, `composition.csv`
(multi-component only), `density_profile.csv`, `summary.txt`.

Interpreting it:

- **Rg vs size** is the compactness check. For a compact globule Rg and N^(1/3) move
  together; if Rg grows faster than N^(1/3) the aggregate is open or elongated.
- **rel_anisotropy** near 0 means spherical; 0.5–0.8 is rod-like, which for a stack of
  planar molecules usually means a one-dimensional slip-stacked column.
- **density profile peaking away from r = 0** is the signature of a hollow structure
  (vesicle); a monotonic decay means a solid core.
- **Composition drifting over time** in a co-assembly indicates demixing — compare the
  first and last quarters of `composition.csv`, not a single frame.

Two guards worth remembering: the cluster criterion is a choice (contact vs centroid,
`--cutoff`), and with a percolating network the "largest cluster" can be most of the box
without being a physical aggregate. State the criterion whenever you quote a size.

### 5b. Stacking geometry — only for planar aromatic systems

```powershell
& $py "$S\analyze.py" --dir <dir> --traj 04_prod.dcd --last-ns 50
& $py "$S\analyze.py" --dir <dir> --traj nowater.dcd --frame-ps 124 --plane-mode intra
& $py "$S\analyze.py" --dir <dir> --traj 04_prod.dcd --outdir <writable path> --max-frames 500
```

Works on a full-system or a solute-only trajectory (the atom count decides which; a
solute-only file must be residue-ordered). Outputs to `<dir>\analysis` by default:
`rdf_<a>_<b>.csv`, `planes_<A>_<B>_angle.csv`, `planes_<A>_<B>_distance.csv`,
`clusters.csv`, `summary.txt` — CSV on purpose, so the numbers go straight into Origin
or Grace.

What it computes:

- **RDF** between two sites, intermolecular pairs only (same-molecule pairs are always
  excluded), per-frame volume and minimum image, `g(r)` normalised by `(N-1)/V`. A site is
  either an atom name (`C5`) or a ring centroid (`ring:ring1`), so the same command works
  on any molecule. `--rdf "C5|O3,ring:ring1|ring:ring1"` overrides the pairs; `A:B` still
  works for plain atom names.
- **Plane geometry** via best-fit SVD normals per molecule, with the planarity rms
  reported so a bad atom group shows up immediately. `--plane-mode inter` (default)
  compares groups on *different* molecules within `--plane-cutoff`; `--plane-mode intra`
  compares the two groups inside one molecule. Interplanar distance is the perpendicular
  component of the centroid separation.
- **Slip angle and J/H classification.** The angle between the two plane *normals* is near
  zero for both slipped and face-to-face stacks, so it cannot tell J from H. What does is
  the **slip angle**: the angle between the centroid-centroid vector and the plane normal
  (0° = face-to-face H-type, larger = displaced J-type). Output
  `planes_<A>_<B>_slip.csv`; the split uses `--slip-threshold` (default 32°).
- **Aggregates** by union-find, `--cluster-mode contact` (true minimum atom distance,
  default) or `centroid` (cheaper, cruder), with `clusters.csv` for the size distribution
  and `clusters_timeseries.csv` for how the assembly evolves frame by frame. The summary
  reports `growing` / `stable` / `shrinking` by comparing the first and last quarter.

The atom groups are whatever `inspect_molecule.py` wrote into `case.json` — there are no
built-in defaults for any particular molecule, deliberately. If `case.json` has no groups
and no `--planes/--rdf` is passed, the script says so instead of guessing, and only the
cluster analysis runs.

### 6. Are the sibling systems comparable?

```powershell
& $py "$S\compare_systems.py" --root <parent dir holding one system per subdirectory>
& $py "$S\compare_systems.py" --root <parent> --resname LIG --json
```

One row per system — molecules, atom count, box, **production-time box volume and
concentration**, production ns, frame spacing — followed by anything that blocks a fair
comparison: concentration spread beyond 5 %, mismatched frame spacing or run length,
differing `1-4scaling`/`timestep`/`cutoff`, differing molecule counts, and (per system)
a knob that changes between its own equilibration and production.

The concentration comes from the **last production frame's unit cell**, not the initial
`CRYST1`: NPT relaxes the box during equilibration, so the starting volume says little
about the conditions the aggregates actually formed under. Water-stripped VMD exports
(`nowater.dcd`, `*_nowater.dcd`, `small.dcd`) are excluded — they are not production runs.

### 7. Can the study support its claims yet?

```powershell
& $py "$S\study_report.py" --root <parent dir> --analysis-root <where analyze.py wrote>
```

One row per system (concentration, production ns, J-type %, median interplanar d,
nearest-neighbour distance, largest cluster, assembly trend) and then a verdict per
claim, drawn from `reference/protocol.md`:

- `J-type slip stacking exists` — every system above 50 % J-type;
- `Interplanar distance <= 4 A` — median of the stacked pairs at or below 4.0 Å;
- `Assembly has equilibrated` — no system still growing, otherwise a size ranking is
  measuring the clock rather than the chemistry;
- `Substituent comparison` — blocked while concentration, frame spacing, run length or
  `oneFourScaling` disagree between systems.

Exit code 2 means at least one claim is blocked. Quoting a cross-system ranking while
that last line is BLOCKED is the mistake this command exists to prevent.

## Reading the results

- A trajectory whose header NSET equals the frames counted on disk is intact up to its
  last frame; `partial_bytes` is only the frame that was being written.
- Two files can both be complete and still not concatenate: compare `last_step + nsavc`
  of one with `first_step` of the next, which is what `inspect_traj.py` does.
- `INTERMOLECULAR` matters: an RDF that includes same-molecule pairs shows a spurious
  peak at the intramolecular distance. Check that the reported peak positions are not
  the intramolecular ones (inside one molecule, the two chosen atoms sit at a fixed short distance).
- A cluster histogram of `50, 1.0` means every frame had one aggregate containing all
  molecules. With a contact cutoff of 4 A that also happens when a percolating network of
  touching molecules links the whole box, which is **not** the same claim as "one
  aggregate" — before reporting cluster sizes, state the criterion (contact vs centroid,
  the cutoff) and check `clusters_timeseries.csv` for whether it changes during the run.
- Self-assembly claims split into two kinds and only one survives a concentration
  difference: stacking **motifs** (J vs H, interplanar distance, which rings pair) are
  local and comparable across concentrations; aggregate **size** and kinetics follow
  concentration directly.
- Before attributing an aggregate-size difference to the substituent, check the
  concentration column of `compare_systems.py`. Two systems can differ by 1.5x in solute
  concentration purely because of how far Packmol spread them, and aggregate size follows
  concentration monotonically.
- A production `.conf` whose `1-4scaling` differs from its own `02_nvt`/`03_equil_npt`
  introduces a potential-energy step at the handoff; the velocities carried over no longer
  correspond to the new Hamiltonian.

## Environment notes

- NAMD is usually not installed locally in this project; scripts here prepare inputs and
  read outputs, they never launch NAMD. Print the command instead of running it:
  `namd2 <new>.conf > <new>.log`.
- Trajectories are large; `analyze.py` streams frame by frame and only reads the atom
  indices it needs, but a 240k-atom DCD at 10 ps/frame over 50 ns is still ~27 GB of
  reads. Use `--max-frames` or `--stride` while iterating.
- If the case directory is read-only for the agent, keep the analysis outputs elsewhere
  with `--outdir` (absolute paths are allowed).
