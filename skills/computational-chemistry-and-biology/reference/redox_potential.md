# Redox potential

## Input

A 3D structure (`StructureInput`): a `rowan.Molecule`, `stjames.Molecule`, or RDKit `Mol` carrying coordinates. Get one any way: load coordinates with `rowan.Molecule.from_xyz_file(path)`, reuse a prior result's `.molecule`, or embed from a SMILES with `rowan.Molecule.from_smiles(...)`.

## Example

```python
import rowan

folder = rowan.get_folder("examples")

wf = rowan.submit_redox_potential_workflow(
    initial_molecule=rowan.Molecule.from_smiles("C1=CC=C(C=C1)C(=O)O"),
    oxidation=True,
    reduction=False,
    folder=folder,
)

result = wf.result()
print(result)   # oxidation and/or reduction potentials, in volts vs SCE
```

## Settings

- `oxidation` (default `True`) and `reduction` (default `False`): which potential(s) to compute. Set both `True` for both.
- `mode` (default `rapid`): the level of theory, trading speed for accuracy across `reckless`, `rapid`, `careful`, and `meticulous`. `rapid` is usually sufficient; step up to `careful` or `meticulous` for higher-accuracy DFT, or `reckless` for a fast, rough estimate.

Predicted potentials are referenced to the saturated calomel electrode (SCE). In `rapid` mode, the workflow reaches a mean absolute error of about 0.32 V on the OROP benchmark set.
