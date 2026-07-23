# Docking

## Input

A protein, a binding pocket, and a single ligand.

- Protein: a `rowan.Protein` or its UUID. Get one from the PDB with `rowan.create_protein_from_pdb_id(pdb_code)`, or upload your own PDB with `rowan.upload_protein(name, path)`. Call `protein.prepare()` to fix nonstandard residues, add missing atoms, and add hydrogens.
- `pocket`: the search box as two `[x, y, z]` points, `[[center_x, center_y, center_z], [size_x, size_y, size_z]]`, the box center and its dimensions in angstroms.
- `initial_molecule`: the ligand as a 3D structure (`StructureInput`): a `rowan.Molecule`, `stjames.Molecule`, or RDKit `Mol` carrying coordinates. Get one any way: embed from a SMILES with `rowan.Molecule.from_smiles(...)`, load coordinates with `rowan.Molecule.from_xyz_file(path)`, or reuse a prior result's `.molecule`. Vina re-poses the ligand in the box, so the input coordinates are only a starting point.

Docking generates poses with AutoDock Vina and refines them with an NNP strain-energy correction by default. Optional conformer search and geometry optimization can run on the ligand beforehand.

## Blind docking (unknown pocket)

`pocket` is always required, but you can dock blind straight from this workflow: set `pocket` to a box large enough to span the whole protein (a large `size` centered on the structure). Use `executable="qvina-w"` for blind docks, since QVina-W is optimized for them; it does not support `vinardo`, so pair it with `scoring_function="vina"`.

To narrow a blind search to likely sites first, run pocket detection and feed a detected pocket straight into `pocket`:

```python
pockets = rowan.submit_pocket_detection_workflow(protein).result().pockets
best = max(pockets, key=lambda p: p.score)  # rank by druggability score
pocket = [list(best.pocket_center), list(best.pocket_sides)]  # ready for submit_docking_workflow
```

## Example

```python
import rowan

folder = rowan.get_folder("examples")

# Dasatinib redocked into its ABL1 co-crystal structure (PDB: 2GQG)
dasatinib = rowan.Molecule.from_smiles("Cc1nc(Nc2ncc(C(=O)Nc3c(C)cccc3Cl)s2)cc(N2CCN(CCO)CC2)n1")

protein = rowan.create_protein_from_pdb_id("2GQG")  # warns if multiple chains
if len(protein.chains) > 1:
    protein = protein.select_chains([protein.chains[0]])
protein.prepare()

center = [44.59, 79.75, 39.59]
size = [24.15, 21.33, 19.88]
wf = rowan.submit_docking_workflow(
    protein,
    pocket=[center, size],
    initial_molecule=dasatinib,
    folder=folder,
)

result = wf.result()
for score in result.scores:
    print(score.score, score.mmgbsa_score, score.posebusters_valid)
```

Each pose's `strain` is its energy above the ligand's lowest-energy conformer, populated only when `do_csearch=True` (otherwise `None`). The result also exposes `best_pose` (the top-scoring pose) and `conformers`. Poses come back with explicit hydrogens reconstructed (Vina strips them during docking), so `best_pose` is a complete 3D structure ready for downstream use such as MD.

Each pose also has an optional `mmgbsa_score`, an MM/GBSA binding free energy estimate in
kcal/mol. Treat it as optional and check for `None`; the result schema does not guarantee that
every docking run or pose has an MM/GBSA value.

## Settings

- `docking_settings`: a `rowan.VinaSettings` or `rowan.GninaSettings` object controlling the docking engine. If provided, it overrides `executable`, `scoring_function`, `exhaustiveness`, and `max_poses` below.
- `executable` (default `vina`): Vina docking implementation, one of `vina`, `qvina2`, `qvina-w`.
- `scoring_function` (default `vinardo`): Vina scoring function, `vinardo` or `vina`. Vinardo is more accurate; Vina is faster. QVina implementations (`qvina2`, `qvina-w`) do not support `vinardo`, so switch to `vina` scoring if you select one.
- `exhaustiveness` (default `8`): how many times Vina attempts to find a pose for each conformer. 8 is typical; 32 is relatively careful.
- `max_poses` (default `4`): maximum number of poses generated per input conformer. The total can exceed this when `do_csearch` is on, since each conformer contributes poses.
- `do_csearch` (default `False`): run an OpenConf conformer search on the input before docking, generating an ensemble of starting poses rather than one arbitrary geometry. This is what enables the per-pose `strain` estimate, but it can significantly increase runtime for large systems.
- `do_optimization` (default `False`): run an AIMNet2 optimization on the input ligand before docking. Skip it if the input is already optimized, to save time.
- `do_pose_refinement` (default `True`): run a constrained AIMNet2 optimization on the output poses (gently relieves clashes without erasing the binding mode).

## gnina (noncovalent and covalent docking)

Use covalent docking for ligands that bind to a specific residue (e.g. cysteine-targeting acrylamides)
instead of binding noncovalently — standard docking can't form or break bonds at all.

Pass a `rowan.GninaSettings` object as `docking_settings` to dock with gnina instead of Vina:

```python
gnina_settings = rowan.GninaSettings(scoring_function="gnina_cnn", exhaustiveness=8, max_poses=4)
wf = rowan.submit_docking_workflow(
    protein,
    pocket=[center, size],
    initial_molecule=dasatinib,
    docking_settings=gnina_settings,
    folder=folder,
)
```

- `scoring_function` (default `gnina_cnn`): `gnina_cnn` rescores poses with gnina's convolutional neural network; `vina` disables the CNN and uses standard Vina scoring.
- `exhaustiveness` (default `8`): how many times gnina attempts to find a pose.
- `max_poses` (default `4`): maximum number of poses generated per input conformer.
- `covalent_ligand_atom_index` / `covalent_protein_atom_index`: 0-based, all-atom indices (including hydrogens) of the reacting ligand and protein atoms. Set both together to run covalent docking; leave both unset (the default) for standard noncovalent docking.
