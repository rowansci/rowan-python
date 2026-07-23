# Hydrogen-bond acceptor/donor strength

## Input

A 3D structure (`StructureInput`): a `rowan.Molecule`, `stjames.Molecule`, or RDKit `Mol` carrying coordinates. Get one any way: load coordinates with `rowan.Molecule.from_xyz_file(path)`, reuse a prior result's `.molecule`, or embed from a SMILES with `rowan.Molecule.from_smiles(...)`. A conformer search and geometry optimization run before predicting, so an embedded SMILES conformer is fine.

## Example

```python
import rowan

folder = rowan.get_folder("examples")

wf = rowan.submit_hydrogen_bond_donor_acceptor_strength_workflow(
    initial_molecule=rowan.Molecule.from_smiles("CC(=O)N(C)C"),  # dimethylacetamide
    folder=folder,
)

result = wf.result()
print(result)  # <HydrogenBondDonorAcceptorStrengthResult acceptors=... donors=...>
```

`submit_hydrogen_bond_basicity_workflow` is an alias that submits the same workflow.

## Settings

- `do_csearch` (default `True`): run a conformer search first. Requires `do_optimization=True` (the search conformers must be optimized before prediction); `do_csearch=True` with `do_optimization=False` raises a `ValueError`.
- `do_optimization` (default `True`): optimize structures before predicting.

The acceptor pKBHX and donor pKα strengths are predicted using neural network potentials and r2SCAN-3c DFT.
