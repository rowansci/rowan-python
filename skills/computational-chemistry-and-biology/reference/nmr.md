# NMR

## Input

A 3D structure (`StructureInput`): a `rowan.Molecule`, `stjames.Molecule`, or RDKit `Mol` carrying coordinates. Get one any way: load coordinates with `rowan.Molecule.from_xyz_file(path)`, reuse a prior result's `.molecule`, or embed from a SMILES with `rowan.Molecule.from_smiles(...)`. The geometry is optimized before predicting (so an embedded SMILES conformer is fine); enable `do_csearch` to also search conformers first.

## Example

```python
import rowan

folder = rowan.get_folder("examples")

wf = rowan.submit_nmr_workflow(
    initial_molecule=rowan.Molecule.from_smiles("O[C@H]1[C@H](C(C)C)CC[C@@H](C)C1"),
    do_csearch=True,  # search conformers first (off by default)
    folder=folder,
)

result = wf.result()
for peak in result.predicted_peaks[1]:  # keyed by atomic number: 1 = hydrogen, 6 = carbon
    print(peak)
```

## Settings

- `solvent` (default `"chloroform"`): solvent for the prediction. Must be one of the 11 NMR-supported solvents — chloroform, tetrahydrofuran, dichloromethane, acetone, acetonitrile, dimethylsulfoxide, methanol, water, benzene, toluene, chlorobenzene — anything else raises a `ValueError`.
- `do_csearch` (default `False`): run a conformer search first. The prediction Boltzmann-averages over the conformers found, so enable it for flexible molecules where averaging meaningfully shifts the result; a rigid molecule can skip it. Requires `do_optimization=True` (the search conformers must be optimized before prediction); `do_csearch=True` with `do_optimization=False` raises a `ValueError`.
- `do_optimization` (default `True`): optimize the structure before predicting.
