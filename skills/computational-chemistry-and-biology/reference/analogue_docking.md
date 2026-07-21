# Analogue docking

## Input

A protein, a reference ligand already bound in its pocket, and a list of analogue SMILES.

- `analogues`: a list of SMILES strings for the analogues to dock.
- `analogue_names`: names parallel to `analogues`, used to set each pose's `Molecule.name` for identification (otherwise poses are keyed only by SMILES).
- `initial_molecule`: the reference bound pose, used as the template the analogues are aligned to. Read it from a file with `rowan.Molecule.from_xyz_file(path)`, or pass a `rowan.Molecule`.
- Protein: a `rowan.Protein` or its UUID. Upload your own PDB with `rowan.upload_protein(name, path)`, or get one from the PDB with `rowan.create_protein_from_pdb_id(name, pdb_code, project_uuid=...)`. Call `protein.prepare()` first to fix nonstandard residues, add missing atoms, and add hydrogens.

This workflow generates conformers of each analogue in poses analogous to the bound reference ligand. Local optimization with the docking scoring function and PoseBusters validation are both opt-in (see Settings).

## Getting a bound reference pose

If you do not have a bound pose yet, dock the reference ligand first (see Docking) and reuse its `best_pose`:

```python
bound_pose = rowan.submit_docking_workflow(...).result().best_pose
# a rowan.Molecule in the protein's coordinate frame; pass it as initial_molecule below
```

This is the usual prep for a relative binding free energy perturbation screen: analogue docking places every analogue into the pocket in the protein's coordinate frame, so the ligands are consistently posed and aligned before you build the `rbfe_graph`. Reach for it when your ligands are not yet posed (or not in the same coordinate space as the protein).

## Example

```python
from pathlib import Path
import rowan

folder = rowan.get_folder("examples")
data_dir = Path("examples/data")

analogues = {
    "analogue-1": "CN(C)CCC[C@@]1(c2ccccc2)OCc2cc(C#N)ccc21",
    "analogue-2": "CN(C)CCC[C@@]1(c2ccc(F)cc2)OCc2c(CC)c(C#N)ccc21",
    "analogue-3": "CN(C)CCC[C@@]1(c2ccc(CCC)cc2)OCc2cc(C#N)ccc21",
}

bound_pose = rowan.Molecule.from_xyz_file(str(data_dir / "citalopram_1iep.xyz"))
protein = rowan.upload_protein("1IEP receptor", data_dir / "1iep_receptorH.pdb")

wf = rowan.submit_analogue_docking_workflow(
    analogues=list(analogues.values()),
    analogue_names=list(analogues.keys()),
    protein=protein,
    initial_molecule=bound_pose,
    folder=folder,
)

result = wf.result()
print(result)   # e.g. <AnalogueDockingResult analogues=3 best=(-8.30, 'CN(C)CCC...')>

# result.analogue_scores: docking scores per analogue (keyed by SMILES).
# result.best_poses: top pose per analogue (keyed by SMILES); each pose's .name is the
# analogue name, so re-key by name for a downstream rbfe_graph:
posed = {p.name: p for p in result.best_poses.values()}
```

Each entry in `result.analogue_scores[smiles]` has the Vina `score` and an optional
`mmgbsa_score`, both in kcal/mol. `mmgbsa_score` is populated only when
`run_local_optimization=True`; otherwise it is `None`.

## Settings

- `scoring_function` (default `vinardo`): scoring function, `vinardo` or `vina`. Vinardo is more accurate; Vina is faster.
- `exhaustiveness` (default `8`): how many times Vina attempts to find a pose for each conformer. 8 is typical; 32 is relatively careful.
- `max_poses` (default `4`): maximum number of poses generated per input conformer.
- `num_conformers_per_analogue` (default `100`): maximum number of conformers generated per analogue. 10-50 is suitable for routine use. 100-1000 is recommended when preparing structures for free energy perturbation calculations.
- `require_posebusters` (default `False`): run PoseBusters validity checks and keep only poses that pass.
- `run_local_optimization` (default `False`): optimize each pose within the binding pocket after generation and compute its MM/GBSA binding free energy estimate. Leaving this off is recommended when consistent template alignment is more important, since local optimization can move poses off the template; `mmgbsa_score` is `None` when it is off.
