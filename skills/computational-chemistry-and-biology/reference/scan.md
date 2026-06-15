# Scan

## Input

A 3D structure (`StructureInput`): a `rowan.Molecule`, `stjames.Molecule`, or RDKit `Mol` carrying coordinates. Get one any way: load coordinates with `rowan.Molecule.from_xyz_file(path)`, reuse a prior result's `.molecule`, or embed from a SMILES with `rowan.Molecule.from_smiles(...)`.

## Example

```python
import rowan

folder = rowan.get_folder("examples")

wf = rowan.submit_scan_workflow(
    initial_molecule=rowan.Molecule.from_smiles("O"),
    scan_settings={
        "type": "angle",      # "bond", "angle", or "dihedral"
        "atoms": [2, 1, 3],   # 1-indexed
        "start": 100,
        "stop": 110,
        "num": 5,             # number of points
    },
    calculation_method="gfn2_xtb",
    calculation_engine="xtb",
    folder=folder,
)

result = wf.result()
for coord, energy in result.get_energies(relative=True):
    print(f"{coord:.1f}  {energy:.2f} kcal/mol")
```

## Settings

- `scan_settings` (required): the coordinate to scan. `type` is `bond`, `angle`, or `dihedral`. `atoms` lists the 1-indexed atoms (2 for a bond, 3 for an angle, 4 for a dihedral). `start` and `stop` give the range, and `num` the number of points. Pass a **list** of these dicts to scan several coordinates *simultaneously* (a concerted scan); every coordinate in the list must use the same `num`.
- `scan_settings_2d` (default none): an optional second coordinate (or list of coordinates), same form as `scan_settings`. Supplying it makes a **2D grid** — every `scan_settings` value is scanned against every `scan_settings_2d` value, giving `num × num_2d` points.
- `calculation_method` (default `omol25_conserving_s`) and `calculation_engine` (default auto-selected): level of theory for each constrained optimization along the scan.
- `mode` (default `auto`): geometry-optimization mode setting the convergence thresholds for each constrained optimization, one of `auto`, `reckless`, `rapid`, `careful`, `meticulous`. `auto` normally resolves to `rapid`. Tighter modes converge each point more strictly at higher cost.
- `constraints` (default none): additional geometric constraints held fixed at every scan point, beyond the scanned coordinate — same as in a basic calculation. Pass `rowan.Constraint(constraint_type="bond", atoms=[1, 2], value=1.5)` (types: `bond`, `angle`, `dihedral`, `freeze_atoms`; atoms 1-indexed; omit `value` to hold the current geometry). Atom order matters for angles and dihedrals (A–B–C is not the same as A–C–B).
- `wavefront_propagation` (default `True`): seed each scan point from already-solved neighboring points. This greatly reduces the influence of the starting geometry and yields a much cleaner surface, at roughly double the runtime. Keep it on unless you have a specific reason not to.

For a 2D scan, `result.get_energies()` returns energies indexed by scan point (the grid is flattened) rather than paired with a single scalar coordinate.

```python
# 2D grid: scan a dihedral against a bond length.
wf = rowan.submit_scan_workflow(
    initial_molecule=mol,
    scan_settings={"type": "dihedral", "atoms": [1, 2, 3, 4], "start": 0, "stop": 180, "num": 10},
    scan_settings_2d={"type": "bond", "atoms": [2, 3], "start": 1.3, "stop": 1.6, "num": 7},
    folder=folder,
)
```

## Result fields

- `get_energies(relative=False)`: list of `(coordinate, energy)` tuples — the scanned coordinate value paired with its energy in Hartree, or in kcal/mol relative to the lowest point when `relative=True`. For a 2D grid the coordinate is the flattened scan-point index.
- `scan_points`: the `Calculation` object at each scan point.
- `scan_point_uuids`: UUIDs of those calculations.
- `messages`: any messages or warnings from the run.
