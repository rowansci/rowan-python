# Interaction energy decomposition

## Input

A single 3D structure (`StructureInput`: a `rowan.Molecule`, `stjames.Molecule`, or RDKit `Mol` with coordinates) holding one non-covalent complex of two interacting fragments, typically read from coordinates with `rowan.Molecule.from_xyz(...)`. It is one system, not two separate molecules: assign one fragment with `fragment1_indices` (1-indexed), and the remaining atoms form the second fragment. The two fragments interact non-covalently, so the split falls between them, not across a covalent bond. SAPT0 then decomposes the interaction energy between the two fragments.

## Example

```python
import rowan

folder = rowan.get_folder("examples")

# A non-covalent dimer loaded from a file; atoms 1-18 are one fragment.
dimer = rowan.Molecule.from_xyz_file("dimer.xyz")

wf = rowan.submit_interaction_energy_decomposition_workflow(
    initial_molecule=dimer,
    fragment1_indices=list(range(1, 19)),  # atoms 1-18 form fragment 1
    folder=folder,
)

result = wf.result()
print(
    result.total_interaction_energy
)  # kcal/mol; also electrostatic/exchange/dispersion/induction components
```

## Settings

- `fragment1_indices` (required): 1-indexed atom indices defining the first fragment. The remaining atoms form the second fragment.
- `method` (default `"sapt0"`): the only supported method. The interaction energy is decomposed into SAPT0 components (electrostatics, exchange, induction, dispersion).
- `basis_set` (default `"jun-cc-pVDZ"`): basis set for the SAPT0 calculation.
