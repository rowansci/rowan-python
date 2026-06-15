# pKa

## Input

The required input type depends on the `method`:

- The two Rowan 3D methods, `aimnet2_wagen2024` and `gxtb_wagen2026`, require a 3D structure with coordinates (a `rowan.Molecule`, `stjames.Molecule`, or RDKit `Mol`), for example `rowan.Molecule.from_smiles("c1ccccc1O")`. Passing a SMILES string raises `ValueError`.
- The two SMILES methods, `chemprop_nevolianis2025` and `starling`, require a SMILES string passed as `initial_molecule`. Passing a 3D structure raises `ValueError`. To use a structure you already have, extract its SMILES first (`molecule.smiles` for a `rowan.Molecule`).

Predicts microscopic (per-site) pKa values: the workflow enumerates protonation and deprotonation microstates and reports a pKa for each.

## Example

```python
import rowan

folder = rowan.get_folder("examples")

wf = rowan.submit_pka_workflow(
    initial_molecule=rowan.Molecule.from_smiles("c1ccccc1O"),
    method="gxtb_wagen2026",
    folder=folder,
)

result = wf.result()
print(result.strongest_acid)     # strongest acidic site pKa, or None
print(result.strongest_base)     # strongest basic site pKa, or None
for ms in result.conjugate_bases:
    print(ms.atom_index, ms.pka)  # per-site deprotonation pKa values
```

## Settings

- `method` (default `gxtb_wagen2026`): algorithm used to predict pKa.
  - `aimnet2_wagen2024`: AIMNet2-based. Requires a 3D structure. Water only.
  - `gxtb_wagen2026`: g-xTB-based. Requires a 3D structure. Water only. Covers the full periodic table.
  - `chemprop_nevolianis2025`: Chemprop-based. Requires SMILES. Supports several solvents. Trained on deprotonation data only, so `protonate_elements` is disabled.
  - `starling`: SMILES-based. Water only.
- `pka_range` (default `(2, 12)`): range of pKa values to report.
- `solvent` (default `"water"`): solvent for the prediction. Only `chemprop_nevolianis2025` accepts non-water solvents (water, DMSO, DMF, acetonitrile, methanol, ethanol, ethylene glycol, N-methylpyrrolidone). The other three methods raise `ValueError` for anything but water.
- `deprotonate_elements` (default `[7, 8, 16]`): atomic numbers of elements whose protons are candidates for deprotonation (nitrogen, oxygen, sulfur).
- `protonate_elements` (default `[7]`): atomic numbers of elements that are candidates for protonation (nitrogen). Not supported with `chemprop_nevolianis2025`, which was trained on deprotonation data only: it is forced empty there, and passing a non-empty value raises a `ValueError`.
- `mode` (default `careful`): thoroughness of the conformer search for the 3D methods (how many conformers are generated and kept, and the energy and RMSD cutoffs). `careful` gives the best balance of speed and accuracy; `rapid` and `reckless` trade accuracy for throughput.

## Result fields

- `strongest_acid`, `strongest_base`: the most acidic and most basic site pKa values, or `None`.
- `conjugate_acids`, `conjugate_bases`: lists of microstates. Each has `atom_index` and `pka`. SMILES-based methods also populate `smiles`. The 3D methods populate `delta_g` (kcal/mol). `chemprop_nevolianis2025` also populates `uncertainty`.
- `structures`: optimized structure calculations, available only for the 3D structure-based methods.

The `aimnet2_wagen2024` method reaches mean absolute errors of around 1 pKa unit on benchmarks such as SAMPL7.
