# Ion mobility

## Input

A 3D structure (`StructureInput`): a `rowan.Molecule`, `stjames.Molecule`, or RDKit `Mol` carrying coordinates. Get one any way: load coordinates with `rowan.Molecule.from_xyz_file(path)`, reuse a prior result's `.molecule`, or embed from a SMILES with `rowan.Molecule.from_smiles(...)`. A conformer search and optimization run first, so an embedded SMILES conformer is fine. The input must be a single connected molecule; mixtures are rejected to avoid ambiguous CCS assignments.

## Example

```python
import rowan

folder = rowan.get_folder("examples")

wf = rowan.submit_ion_mobility_workflow(
    initial_molecule=rowan.Molecule.from_smiles("c1ccccn1"),
    protonate=True,  # model the protonated ion that positive-mode IM-MS detects
    folder=folder,
)

result = wf.result()
print(result.average_ccs)  # weighted-mean CCS in A^2 (also average_ccs_stdev for the uncertainty)
```

## Settings

- `temperature` (default `300`): temperature in K.
- `protonate` (default `False`): enumerate protomers by protonating the N, O, and S sites, then run the CCS calculation on the lowest-energy protomer. This models the protonated ion (the molecule plus one proton, written [M+H]+ in mass spectrometry, where M is the molecule) that positive-mode ion-mobility MS detects.
- `do_csearch` (default `True`): run a conformer search first. Requires `do_optimization=True` (the search conformers must be optimized before the CCS calculation); `do_csearch=True` with `do_optimization=False` raises a `ValueError`.
- `do_optimization` (default `True`): optimize before computing CCS.
