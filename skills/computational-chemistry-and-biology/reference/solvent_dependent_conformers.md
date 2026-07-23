# Solvent-dependent conformers

## Input

A molecule with 3D coordinates. Generate one from a SMILES with `rowan.Molecule.from_smiles(...)` (RDKit embeds it into 3D), or pass existing coordinates with `rowan.Molecule.from_xyz(...)`.

## Example

```python
import rowan

folder = rowan.get_folder("examples")

wf = rowan.submit_solvent_dependent_conformers_workflow(
    initial_molecule=rowan.Molecule.from_smiles("CC(=O)N[C@@H](C)C(=O)NC"),  # alanine dipeptide
    folder=folder,
)

result = wf.result()
print(result)  # <SolventDependentConformersResult conformers=... solvents=...>

# Transfer free energies relative to the lowest solvent
for solvent, dg in result.relative_free_energy_by_solvent.items():
    print(f"  {solvent.value}: dG = {dg:.2f} kcal/mol")

# Per-conformer populations in water
for i, conf in enumerate(result.conformers):
    pop = conf.population_by_solvent.get(rowan.Solvent.WATER, 0)
    print(f"  Conformer {i + 1}: population={pop * 100:.1f}% (water)")
```

## Settings

- `solvents` (default hexane, octanol, chloroform, DMSO, and water): solvents to score conformers in, as a list of `rowan.Solvent` enum values. Scoring uses CPCM-X.
- `conf_gen_settings`: conformer generation settings. When omitted, inherits the stjames default for this workflow (currently OpenConf). Other options: `rowan.ETKDGSettings`, `rowan.iMTDSettings`, `rowan.iMTDGCSettings`. The energy window lives on this object (e.g. `OpenConfSettings.energy_window_kcal`, default 10), not as a separate argument.

## Result fields

- `relative_free_energy_by_solvent`: transfer free energy per solvent, relative to the lowest (kcal/mol).
- `conformers`: per-conformer results, each with `population_by_solvent` and `relative_free_energy_by_solvent`.
- `per_solvent_properties`: per-solvent ensemble averages — `solvent_accessible_surface_area`, `polar_solvent_accessible_surface_area`, `radius_of_gyration`.
