# Descriptors

## Input

A 3D structure (`StructureInput`): a `rowan.Molecule`, `stjames.Molecule`, or RDKit `Mol` carrying coordinates. Get one any way: load coordinates with `rowan.Molecule.from_xyz_file(path)`, reuse a prior result's `.molecule`, or embed from a SMILES with `rowan.Molecule.from_smiles(...)`. The geometry is optimized first, so an embedded SMILES conformer is fine.

## Example

```python
import rowan

folder = rowan.get_folder("examples")

wf = rowan.submit_descriptors_workflow(
    initial_molecule=rowan.Molecule.from_smiles("CC(=O)Oc1ccccc1C(=O)O"),  # aspirin
    solvent="water",   # opt into the COSMO descriptors; omit for gas phase
    folder=folder,
)

result = wf.result()
print(result)               # <DescriptorsResult n=...>
print(result.descriptors)   # dict of computed descriptors
```

`result.descriptors` holds electronic descriptors from xTB (atomic charges, the global electrophilicity index, and Fukui indices, as in the Fukui workflow), cheminformatic descriptors from Mordred, and the COSMO solvation descriptors when a `solvent` is set.

## Settings

- `solvent` (default none, gas phase): solvent used for the COSMO descriptors (surface area, screening charge, dielectric energy, polar surface area). When provided, the additional COSMO descriptors are computed; left at none, the workflow runs in the gas phase and skips them. Introspect `rowan.Solvent` for the available solvents.
- `do_optimization` (default `True`): run a GFN2-xTB geometry optimization before computing descriptors. Set it `False` only when your input is already optimized at a comparable level, to skip redundant work; the descriptors are geometry-dependent, so a poor geometry gives poor descriptors.
