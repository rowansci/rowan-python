# Spin states

## Input

A 3D structure (`StructureInput`): a `rowan.Molecule`, `stjames.Molecule`, or RDKit `Mol` carrying coordinates. Get one any way: load coordinates with `rowan.Molecule.from_xyz_file(path)`, reuse a prior result's `.molecule`, or embed from a SMILES with `rowan.Molecule.from_smiles(...)`. The same optimization settings are applied to each spin state, and the lowest-energy one is reported.

## Example

```python
import rowan

folder = rowan.get_folder("examples")

wf = rowan.submit_spin_states_workflow(
    initial_molecule=rowan.Molecule.from_smiles("[CH2]"),
    states=[1, 3],  # singlet vs triplet; defaults to [1, 3, 5] if omitted
    folder=folder,
)

result = wf.result()
print(result)  # e.g. <SpinStatesResult states=2 ground=(mult=3, E=-39.14 H)>

# Energies relative to the ground (lowest-energy) state.
for state, rel in zip(result.spin_states, result.get_energies(relative=True)):
    print(f"  multiplicity {state.multiplicity}: {rel:.2f} kcal/mol")
```

## Settings

- `states` (default `[1, 3, 5]`, or `[2, 4, 6]` for odd-electron molecules): list of spin multiplicities to compare, e.g. `[1, 3]` for singlet vs triplet. All states must share the molecule's electron parity.
- `multistage_opt_settings` (default none): a `MultiStageOptSettings` applied to every state, bundling one or more `optimization_settings` stages and a `singlepoint_settings`. Omit to use the default `r2scan_3c//gfn2_xtb` stack. To run in implicit solvent, set `solvent_settings` on the `singlepoint_settings` — optimizations stay gas-phase and the solvent model is applied only to the final energy.
- `frequencies` (default `False`): compute vibrational frequencies on the final optimization of each state.
- `transition_state` (default `False`): optimize each state to a transition state rather than a minimum.
- `constraints` (default none): geometric constraints (`bond`, `angle`, `dihedral`, or `freeze_atoms`) held fixed during every state's optimization, as in a scan or basic calculation. Pass `rowan.Constraint(constraint_type="bond", atoms=[1, 2], value=1.5)`.

## Result fields

- `spin_states`: states in submission order, each with `multiplicity`, `energy` (Hartree), and `calculation_uuids`.
- `get_energies(relative=False)`: per-state energies in Hartree, or kcal/mol relative to the ground state when `relative=True`.
- `get_calculation(multiplicity, stage=-1)`: fetch the `Calculation` for one state (default final stage).
- `messages`: any messages or warnings from the run.
