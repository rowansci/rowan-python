# Docking

## Input

A protein, a binding pocket, and a single ligand.

- Protein: any stored `rowan.Protein` or protein UUID. Get one from the PDB with `rowan.create_protein_from_pdb_id(pdb_code)`, or upload your own PDB with `rowan.upload_protein(name, path)`. Protein preparation is recommended before docking; when chaining from it, passing `prepared_protein_uuid` avoids fetching structure data solely for submission.
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
protein = protein.select_chains(["A"])
preparation_workflow = rowan.submit_protein_preparation_workflow(
    protein=protein.uuid, folder=folder
)
prepared_protein_uuid = preparation_workflow.result().prepared_protein_uuid

center = [44.59, 79.75, 39.59]
size = [24.15, 21.33, 19.88]
wf = rowan.submit_docking_workflow(
    prepared_protein_uuid,
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

## gnina docking

`GninaSettings` selects gnina, but does not by itself make the run covalent. With neither covalent
atom index set, gnina performs standard noncovalent docking. Set both indices to form a bond between
a known ligand atom and protein atom.

For noncovalent gnina docking with CNN scoring:

```python
gnina_settings = rowan.GninaSettings(scoring_function="gnina_cnn", exhaustiveness=8, max_poses=4)
wf = rowan.submit_docking_workflow(
    protein.uuid,
    pocket=[center, size],
    initial_molecule=dasatinib,
    docking_settings=gnina_settings,
    folder=folder,
)
```

Use covalent docking for a ligand designed to bind a known residue, such as a cysteine-targeting
acrylamide. Prepare the protein first and find the reactive protein atom in the prepared structure,
because preparation can change atom ordering and residue numbering:

Supply the ligand in its expected post-reaction, covalently bound topology; gnina does not infer the
reaction. For a Michael acceptor `C=CC(=O)NR`, use the hydrogen-capped product `CCC(=O)NR` and
select the terminal β-carbon as the covalent ligand atom.

```python
prepared_protein = preparation_workflow.result().get_prepared_protein()
reactive_protein_atom_index = prepared_protein.get_atom_index(
    chain="A", residue=reactive_residue, atom="SG"
)
gnina_settings = rowan.GninaSettings(
    scoring_function="vina",
    covalent_ligand_atom_index=reactive_ligand_atom_index,
    covalent_protein_atom_index=reactive_protein_atom_index,
)
wf = rowan.submit_docking_workflow(
    prepared_protein.uuid,
    pocket=[center, size],
    initial_molecule=ligand,
    docking_settings=gnina_settings,
    folder=folder,
)
```

Both atom indices are zero-based and include hydrogens. The ligand index refers to the atom order in
`initial_molecule`; the protein index refers to PDB atom-record order. For XYZ and SDF inputs, atom
indices follow file order. For `Molecule.from_smiles`, heavy atoms follow SMILES parse order
(`CCC(=O)N`: C(0)-C(1)-C(2)-O(3)-N(4)), with hydrogens appended afterward, so the reactive
heavy-atom index can normally be counted directly. Only unusual SMILES containing explicit `[H]`
atoms require verifying the index with RDKit's `GetIdx()`; bracket hydrogen counts such as `[C@H]`
do not create separate atoms.

Covalent mode requires `scoring_function="vina"`; `gnina_cnn` is available only for noncovalent
gnina docking.
PoseBusters validation is skipped for covalent poses, so `posebusters_valid` is `None` and means not
evaluated rather than failed. Do not use it to reject covalent poses.

- `scoring_function` (default `gnina_cnn`): `gnina_cnn` rescores poses with gnina's convolutional neural network; `vina` disables the CNN and uses standard Vina scoring.
- `exhaustiveness` (default `8`): how many times gnina attempts to find a pose.
- `max_poses` (default `4`): maximum number of poses generated per input conformer.
- `covalent_ligand_atom_index` / `covalent_protein_atom_index`: 0-based, all-atom indices
  (including hydrogens) of the reacting ligand and protein atoms. Set both together to run covalent
  docking; leave both unset for noncovalent docking.
