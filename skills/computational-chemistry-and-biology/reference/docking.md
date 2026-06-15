# Docking

## Input

A protein, a binding pocket, and a single ligand.

- Protein: a `rowan.Protein` or its UUID. Get one from the PDB with `rowan.create_protein_from_pdb_id(name, pdb_code, project_uuid=...)`, or upload your own PDB with `rowan.upload_protein(name, path)`. Call `protein.prepare()` first to fix nonstandard residues, add missing atoms, and add hydrogens.
- `pocket`: the search box as two `[x, y, z]` points, `[[center_x, center_y, center_z], [size_x, size_y, size_z]]`, the box center and its dimensions in angstroms.
- `initial_molecule`: the ligand as a 3D structure (`StructureInput`): a `rowan.Molecule`, `stjames.Molecule`, or RDKit `Mol` carrying coordinates. Get one any way: embed from a SMILES with `rowan.Molecule.from_smiles(...)`, load coordinates with `rowan.Molecule.from_xyz_file(path)`, or reuse a prior result's `.molecule`. Vina re-poses the ligand in the box, so the input coordinates are only a starting point.

Docking generates poses with AutoDock Vina and refines them with an NNP strain-energy correction by default. Optional conformer search and geometry optimization can run on the ligand beforehand.

## Blind docking (unknown pocket)

`pocket` is always required, but you can dock blind straight from this workflow: set `pocket` to a box large enough to span the whole protein (a large `size` centered on the structure). Use `executable="qvina-w"` for blind docks, since QVina-W is optimized for them; it does not support `vinardo`, so pair it with `scoring_function="vina"`.

To narrow a blind search to likely sites first, run pocket detection and feed a detected pocket straight into `pocket`:

```python
pockets = rowan.submit_pocket_detection_workflow(protein).result().pockets
best = max(pockets, key=lambda p: p.score)   # rank by druggability score
pocket = [list(best.pocket_center), list(best.pocket_sides)]   # ready for submit_docking_workflow
```

## Example

```python
import rowan

folder = rowan.get_folder("examples")

protein = rowan.create_protein_from_pdb_id(
    "CDK2", "1HCK", project_uuid=rowan.default_project().uuid
)
protein.prepare()

wf = rowan.submit_docking_workflow(
    protein.uuid,
    pocket=[[103.55, 100.59, 82.99], [27.76, 32.67, 48.79]],  # [center], [box size] in Angstrom
    initial_molecule=rowan.Molecule.from_smiles("CCC(C)(C)NC1=NCC2(CCC(=O)C2C)N1"),
    folder=folder,
)

result = wf.result()
for score in sorted(result.scores, key=lambda s: s.score):
    print(score.score, score.posebusters_valid, score.strain)   # score, PoseBusters validity, strain
```

Each pose's `strain` is its energy above the ligand's lowest-energy conformer, populated only when `do_csearch=True` (otherwise `None`). The result also exposes `best_pose` (the top-scoring pose) and `conformers`. Poses come back with explicit hydrogens reconstructed (Vina strips them during docking), so `best_pose` is a complete 3D structure ready for downstream use such as MD.

## Settings

- `executable` (default `vina`): docking implementation, one of `vina`, `qvina2`, `qvina-w`.
- `scoring_function` (default `vinardo`): scoring function, `vinardo` or `vina`. Vinardo is more accurate; Vina is faster. QVina implementations (`qvina2`, `qvina-w`) do not support `vinardo`, so switch to `vina` scoring if you select one.
- `exhaustiveness` (default `8`): how many times Vina attempts to find a pose for each conformer. 8 is typical; 32 is relatively careful.
- `max_poses` (default `4`): maximum number of poses generated per input conformer. The total can exceed this when `do_csearch` is on, since each conformer contributes poses.
- `do_csearch` (default `False`): run an OpenConf conformer search on the input before docking, generating an ensemble of starting poses rather than one arbitrary geometry. This is what enables the per-pose `strain` estimate, but it can significantly increase runtime for large systems.
- `do_optimization` (default `False`): run an AIMNet2 optimization on the input ligand before docking. Skip it if the input is already optimized, to save time.
- `do_pose_refinement` (default `True`): run a constrained AIMNet2 optimization on the output poses (gently relieves clashes without erasing the binding mode).
