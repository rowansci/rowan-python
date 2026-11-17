# Covalent Inhibitor Scan

> **Beta:** This workflow is in beta. Its interface and behavior may change.

Submit a prepared protein–ligand complex, its ligand reactant SMILES, and required
`rowan.CovalentInhibitorScanSettings`. The prepared structure must contain explicit
hydrogens and the intended protonation and connectivity. Use residue identifiers and
atom names from that structure to select the ligand, reactive region, and coordinates.

```python
import rowan

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
region = [
    rowan.ResidueRef("polymer", "A", f"A.{number}")
    for number in (18, 19, 20, 23, 24, 25, 26, 27, 160, 161, 162, 163, 164, 182, 183, 184, 185)
]
settings = rowan.CovalentInhibitorScanSettings(
    ligand=ligand, reactive_region=[*region, ligand], scan_coordinates=[bond]
)
protein = rowan.upload_protein("Prepared cathepsin K/S2", "prepared-cathepsin-k.pdb")
workflow = rowan.submit_covalent_inhibitor_scan_workflow(
    protein=protein,
    reactant_smiles="Cn1cnc2cnc(C#N)nc21",
    settings=settings,
    is_draft=True,
)
```

These selections describe the prepared cathepsin K/S2 example; adapt them to your
structure. Preparation and choosing a chemically appropriate reactive region remain
caller responsibilities. Call `workflow.submit_draft()` and then `workflow.result()`
to run and retrieve results. The executable `examples/covalent_inhibitor_scan.py`
also demonstrates a second coordinate with `rowan.LinearCombinationCoordinate`.

`result.windows` follows the submitted grid order. Each window contains its grid
index, centers, mean coordinates, sample count, optional reactive-region molecule,
and optional full-system DCD reference. Molecules are shared `rowan.StJamesMolecule`
models; use `rowan.Molecule.from_stjames` for SDK molecule helpers.

`result.surface` contains bin edges, flattened free energies, effective sample counts,
uncertainty, and outside probability. Its axes are ascending, with the last axis
varying fastest, even when the submitted grid descends. Free energies and uncertainties
are in kcal/mol. Missing cells remain `None`; the whole surface can be `None` when
analysis is unavailable. Barrier extraction and convergence assessment are user analysis.
Additional MD restraints remain in the reported surface; only umbrella biases are removed.

Use `result.resolved_setup`, `result.surface_diagnostics`, and `result.messages` to
inspect atom mappings and sampling diagnostics. `result.data` retains the raw result;
use `workflow.fetch_latest().data` for failed or incompatible workflows.

Retrieve the shared solvated topology with `result.get_topology()`. Download its
structure with `topology.download_structure(workflow_uuid=result.workflow_uuid)`.
`result.download_trajectory_files([0, 2], path="trajectories")` downloads selected
DCDs using **positions in `result.windows`**, not the windows' multidimensional
`index` values. Requested order is preserved; selecting a window without a saved DCD
raises an error. Trajectories require API support for covalent window downloads.

::: rowan.workflows.covalent_inhibitor_scan
    handler: python
    options:
      show_source: false
      show_root_heading: false
      show_root_toc_entry: false
      members_order: source
      group_by_category: true
      filters: ["!^_"]
