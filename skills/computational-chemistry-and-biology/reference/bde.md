# BDE

## Input

A 3D structure (`StructureInput`): a `rowan.Molecule`, `stjames.Molecule`, or RDKit `Mol` carrying coordinates. Get one any way: load coordinates with `rowan.Molecule.from_xyz_file(path)`, reuse a prior result's `.molecule`, or embed from a SMILES with `rowan.Molecule.from_smiles(...)`.

## Example

```python
import rowan

folder = rowan.get_folder("examples")

mol = rowan.Molecule.from_smiles("CCO")  # ethanol

# Find the bond you want. The bond-finder helpers return 1-indexed (heavy, H) pairs.
# Here: oxygen (atomic number 8) to hydrogen (1), within 1.2 Angstroms.
oh_bonds = rowan.find_bonds(mol, element_a=8, element_b=1, distance_max=1.2)  # -> [(3, 9)]
h = oh_bonds[0][1]                            # the hydroxyl hydrogen (atom 9)

wf = rowan.submit_bde_workflow(
    initial_molecule=mol,
    fragment_indices=[[h]],   # break the O-H: the dissociating fragment is just that H
    folder=folder,
)

result = wf.result()
print(result.bdes)   # list of BDEEntry, each with fragment_idxs and energy (Hartree)
```

Reported BDEs are empirically corrected to experiment (the ExpBDE54 benchmark) to account for zero-point, enthalpy, and relativistic effects, so they are corrected enthalpies rather than raw electronic energies. Expect a few kcal/mol RMSE versus experiment.

## Settings

- `mode` (default `omol25_conserving_s`): the level of theory, given as a BDE method string: `omol25_conserving_s` (an NNP, the default, recommended for most work), `g_xtb//gfn2_xtb` (semiempirical), or `r2scan3c//gfn2_xtb` (DFT single point on a semiempirical geometry).
- `multistage_opt_settings` (optional): replace the method sequence that `mode` would build, as a `MultiStageOptSettings` (optimization stage(s) plus a final singlepoint). Prefer a `mode` string. Only reach for this when you need a sequence the modes cannot express, such as a singlepoint above r2SCAN-3c.
- Which bonds to break (pick an approach):
  - `fragment_indices=[[...]]`: each inner list is the 1-indexed atoms of one dissociating fragment, which must connect to the rest of the molecule by a single bond. A lone hydrogen `[[9]]` breaks that C-H or O-H; a whole ring breaks a biphenyl-type bond. List several inner lists for several BDEs. To get indices, use the bond finders `rowan.find_bonds(mol, element_a, element_b, distance_max)`, `rowan.find_ch_bonds(mol)`, or `rowan.find_cx_bonds(mol)`, which return 1-indexed atom pairs. If the fragment connects by more than one bond (an internal atom, a ring atom), there is no single bond to break and the BDE comes back `None`.
  - `all_CH=True`: every C-H bond.
  - `all_CX=True`: every C-X bond (carbon to a halogen: F, Cl, Br, I).
