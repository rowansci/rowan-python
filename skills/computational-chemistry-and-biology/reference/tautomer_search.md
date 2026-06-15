# Tautomer search

## Input

A 3D structure (`StructureInput`): a `rowan.Molecule`, `stjames.Molecule`, or RDKit `Mol` carrying coordinates. Get one any way: load coordinates with `rowan.Molecule.from_xyz_file(path)`, reuse a prior result's `.molecule`, or embed from a SMILES with `rowan.Molecule.from_smiles(...)`. Tautomers are enumerated from it and ranked by relative stability.

## Example

```python
import rowan

folder = rowan.get_folder("examples")

wf = rowan.submit_tautomer_search_workflow(
    initial_molecule=rowan.Molecule.from_smiles("O=C1C=CC=CN1"),
    folder=folder,
)

result = wf.result()

best = result.best_tautomer   # predominant tautomer (highest Boltzmann weight), a rowan.Molecule

# All tautomers, with populations and relative energies.
for t in result.tautomers:
    print(f"  weight={t.weight:.2f}, dE={t.predicted_relative_energy:.2f} kcal/mol")
```

## Settings

- `conf_gen_settings` and `multistage_opt_settings`: override conformer generation and the optimization stack. When omitted, both inherit the stjames defaults for this workflow — currently OpenConf (`max_confs=20`) for conformers and an AIMNet2/wB97M-D3 optimization with an AIMNet2/wB97M-D3 CPCMx(water) singlepoint.

## Result fields

- `best_tautomer`: the predominant tautomer (highest Boltzmann weight) as a `rowan.Molecule`.
- `tautomers`: all tautomers, each with `weight` (Boltzmann population), `predicted_relative_energy` (kcal/mol), `energy` (Hartree), and `structure_uuids`.
- `molecules`: the `rowan.Molecule` for every tautomer (one API call per tautomer).
- `messages`: any messages or warnings from the run.

## Accuracy

On the aqueous subset of the TautoBase benchmark, the workflow predicts the correct lowest-energy tautomer about 89% of the time, with a mean absolute error of 2.1 kcal/mol (RMSE 3.0) on relative tautomer energies — comparable to high-level DFT.
