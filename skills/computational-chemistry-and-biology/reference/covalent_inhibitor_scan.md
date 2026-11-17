# Covalent inhibitor scan

> **Beta:** This workflow is under active development. Its interface and behavior may change.

## Input

Pass a prepared `rowan.Protein` or UUID, `reactant_smiles` matching the ligand's
intended protonation and connectivity, and required `rowan.CovalentInhibitorScanSettings`.
The prepared structure must include explicit hydrogens. Preparation and reactive-region
selection remain caller responsibilities.

Settings select the non-polymer `ligand`, whole-residue `reactive_region`, and one or
two `scan_coordinates`. Identify atoms with `rowan.ResidueRef` and
`rowan.ResidueAtomRef` using the prepared structure's identifiers and atom names.
The region must include the ligand and every scan atom. Each `rowan.ScanSettings`
contains a `rowan.DistanceCoordinate` or `rowan.LinearCombinationCoordinate`, plus
`start`, `stop`, `num`, and an optional `force_constant` in kcal/mol/Å².

## Example

```python
import rowan

# These references must match the prepared input structure.
cysteine = rowan.ResidueRef("polymer", "A", "A.25")
ligand = rowan.ResidueRef("non_polymer", "A.300")
bond = rowan.ScanSettings(
    coordinate=rowan.DistanceCoordinate(
        atoms=(rowan.ResidueAtomRef(cysteine, "SG"), rowan.ResidueAtomRef(ligand, "C27"))
    ),
    start=1.6,
    stop=4.6,
    num=31,
    force_constant=100,
)
settings = rowan.CovalentInhibitorScanSettings(
    ligand=ligand,
    reactive_region=selected_region,  # whole residues, including the ligand and scan atoms
    scan_coordinates=[bond],
)
workflow = rowan.submit_covalent_inhibitor_scan_workflow(
    protein=protein,
    reactant_smiles=ligand_smiles,
    settings=settings,
    is_draft=True,
)
```

Review the draft before submitting it. The executable SDK example
`examples/covalent_inhibitor_scan.py` includes a complete region selection for prepared
cathepsin K/S2 and an optional proton-transfer coordinate.

## Settings and results

Simulation and force-field settings live in `CovalentInhibitorScanSettings`, including
`calc_settings`, `equilibration_time_ns`, `simulation_time_ns`, `timestep_fs`,
`frame_save_interval_ps`, and optional `bootstrap_samples`. Additional MD restraints
remain in the free-energy surface. Implicit solvent, pressure control, and REST2
are unsupported.

After submission, `workflow.result()` provides:

- `settings` and `resolved_setup`: submitted settings, atom mappings, capped bonds, charge, and spin.
- `windows`: target centers, mean coordinates, sample counts, optional shared-model molecules, and optional DCD references.
- `surface`: ascending bin edges, flattened free energies in kcal/mol, effective counts, uncertainty, and outside probability. The last axis varies fastest. Missing cells stay `None`; the surface may be unavailable.
- `surface_diagnostics` and `messages`: exchange counts, window overlap, and workflow messages.
- `get_topology()`: shared solvated topology. Download with `workflow_uuid=result.workflow_uuid`.
- `download_trajectory_files(window_indices, path=None)`: DCD files in requested order, selected by positions in `windows`. Missing DCDs raise an error; API support for covalent window downloads is required.

Window order follows the submitted grid and can differ from ascending surface order.
Barrier extraction and convergence assessment are user analysis. Raw data is available
through `result.data`; use `workflow.fetch_latest().data` for failed or incompatible results.
