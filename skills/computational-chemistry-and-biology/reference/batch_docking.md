# Batch docking

## Input

A protein, a binding pocket, and a list of ligand SMILES.

- `smiles_list`: a list of SMILES strings to dock.
- Protein: any stored `rowan.Protein` or protein UUID. Get one from the PDB with `rowan.create_protein_from_pdb_id(pdb_code, name=..., project_uuid=...)`, or upload your own PDB with `rowan.upload_protein(name, path)`. Protein preparation is recommended before docking; when chaining from it, passing `prepared_protein_uuid` avoids fetching structure data solely for submission.
- `pocket`: the search box as two `[x, y, z]` points, `[[center_x, center_y, center_z], [size_x, size_y, size_z]]`, the box center and its dimensions in angstroms.

Each SMILES is prepared internally: RDKit generates 3D coordinates, adds hydrogens, embeds several conformers, and runs a quick MMFF94 optimization, then the lowest-energy conformer is docked. This is why the input is SMILES rather than a posed 3D ligand.

## Example

```python
import rowan

folder = rowan.get_folder("examples")

ligands = [
    "CCC(C)(C)NC1=NCC2(CCC(=O)C2C)N1",
    "CCC(C)CN=C1NCC2(CCCOC2)CN1",
    "CC(C)CCNC1=NCC2CC(COC2=N)O1",
]

protein = rowan.create_protein_from_pdb_id("1HCK", name="CDK2")
preparation_workflow = rowan.submit_protein_preparation_workflow(
    protein=protein.uuid,
    folder=folder,
)
prepared_protein_uuid = preparation_workflow.result().prepared_protein_uuid

wf = rowan.submit_batch_docking_workflow(
    ligands,
    prepared_protein_uuid,
    pocket=[[103.55, 100.59, 82.99], [27.76, 32.67, 48.79]],  # [center], [box size] in Angstrom
    folder=folder,
)

result = wf.result()
print(result.scores)  # dict of SMILES -> best docking score
```

`result.scores` always maps every input SMILES to its best docking score, or `None` when that
ligand failed. `result.refined_scores` provides the same mapping to `rowan.DockingScore`
objects for the top `num_poses_to_save` compounds and `None` for the others. Each saved score
contains the raw pose; with `run_mmgbsa=True`, a successful refinement replaces it with the
refined pose and populates `mmgbsa_score`.

## Settings

- `executable` (default `vina`): docking implementation, one of `vina`, `qvina2`, `qvina-w`.
- `scoring_function` (default `vinardo`): scoring function, `vinardo` or `vina`. Vinardo is more accurate; Vina is faster. The QVina implementations (`qvina2`, `qvina-w`) do not support `vinardo`, so switch to `vina` scoring if you select one.
- `exhaustiveness` (default `8`): how many times the docking engine attempts to find a pose for each ligand. 8 is typical; 32 is relatively careful.
- `num_poses_to_save` (default `0`): number of top-scoring compounds whose best pose to retain. Leave this at zero when only scores are needed; increase it when downstream pose inspection is required.
- `run_mmgbsa` (default `False`): refine the saved poses and compute MM/GBSA binding free energy estimates. This is ignored when `num_poses_to_save=0` and adds cost when enabled, so use it for post-docking reranking rather than routine high-throughput screening.
