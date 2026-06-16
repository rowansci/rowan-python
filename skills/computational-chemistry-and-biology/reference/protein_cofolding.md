# Protein co-folding

## Input

One or more biomolecule sequences, optionally with ligands, supplied directly as sequences and SMILES (no prepared `Protein` object needed). At least one protein, DNA, or RNA sequence is required; ligands alone are not enough, and submitting without a biomolecule sequence raises a `ValueError`.

- `initial_protein_sequences`: a list of protein amino-acid sequences.
- `initial_dna_sequences`, `initial_rna_sequences`: lists of nucleic-acid sequences.
- `initial_smiles_list`: a list of ligand SMILES, all co-folded together with the biomolecules in a single complex. One workflow is one prediction of the whole complex, not one job per ligand; to screen several ligands separately, submit one workflow per ligand.
- `ligand_binding_affinity_index`: index into `initial_smiles_list` of the ligand to predict a binding affinity for.

This workflow predicts 3D structures of biomolecules and protein-ligand complexes using AlphaFold 3-style models (Boltz-2, Boltz-2.1, Chai-1r, Boltz-1, OpenFold-3).

## Example

```python
import rowan

folder = rowan.get_folder("examples")

protein_sequence = "..."  # your target's full amino-acid sequence

wf = rowan.submit_protein_cofolding_workflow(
    initial_protein_sequences=[protein_sequence],
    initial_smiles_list=["CCOC(=O)N1c2ccc(C(F)(F)F)cc2[C@@H](N(Cc2cc(C(F)(F)F)cc(C(F)(F)F)c2)C(=O)OC)C[C@H]1CC"],
    ligand_binding_affinity_index=0,   # predict affinity for the first ligand
    model="boltz_2",
    use_msa_server=True,
    use_potentials=False,
    num_samples=None,
    compute_strain=False,
    do_pose_refinement=True,           # optimize the predicted ligand pose
    folder=folder,
)

result = wf.result()
print(result)   # e.g. <ProteinCofoldingResult predictions=5 iptm=0.87>
```

## Settings

- `model` (default `boltz_2`): co-folding model, one of `boltz_2`, `boltz_2_1`, `chai_1r`, `boltz_1`, `openfold_3`. The model changes which `affinity_score` fields are returned (see Result fields). `boltz_2_1` runs via Boltz's hosted API and is noticeably slower than the locally-run models.
- `use_msa_server` (default `True`): generate a multiple-sequence alignment, which improves co-folding accuracy. Queries run on a secure Rowan-hosted server.
- `use_potentials` (default `False`): use physics-based potentials to steer predictions toward more physical poses (this is what distinguishes Boltz-2x from Boltz-2). Required when a constraint sets `force=True`.
- `num_samples` (default model-dependent): number of predicted samples to generate.
- `compute_strain` (default `False`): run an OpenConf conformer search to estimate the predicted ligand's strain energy. Adds a few minutes (more for large ligands) and requires `do_pose_refinement=True` (`compute_strain=True` with it off raises a `ValueError`).
- `do_pose_refinement` (default `False`): run a constrained AIMNet2 optimization on the output structure to refine the predicted pose.
- `ligand_binding_affinity_index` (default none): index into `initial_smiles_list` of the ligand to predict a binding affinity for.
- `templates` (default none): structural templates that guide the prediction toward a reference protein structure, each a `rowan.CofoldingTemplate(protein=..., max_distance=...)` where `protein` is a PDB or an uploaded protein UUID. Set `max_distance` (in angstroms) to enforce a maximum deviation from the template via potentials (requires `use_potentials=True`). Boltz-2, Boltz-2.1, or OpenFold-3 only.

## Constraints

Encode prior knowledge of where a ligand binds with `pocket_constraints` and `contact_constraints`. Targets are addressed with `rowan.ConstraintTarget(input_type, input_index, token_index)`, where `input_type` is `"protein"` or `"ligand"`, `input_index` selects which entry in the input lists, and `token_index` is the residue index for proteins or the atom index for ligands.

```python
his57 = rowan.ConstraintTarget(input_type="protein", input_index=0, token_index=39)
asp189 = rowan.ConstraintTarget(input_type="protein", input_index=0, token_index=134)
amidine_c = rowan.ConstraintTarget(input_type="ligand", input_index=0, token_index=1)

pocket = rowan.PocketConstraint(
    input_type="ligand",   # the binder being placed into the pocket
    input_index=0,
    contacts=[his57, asp189],
    max_distance=6.0,
)
salt_bridge = rowan.ContactConstraint(
    token_1=amidine_c,
    token_2=asp189,
    max_distance=4.0,
    force=True,            # requires use_potentials=True
)

wf = rowan.submit_protein_cofolding_workflow(
    initial_protein_sequences=[TRYPSIN],
    initial_smiles_list=[BENZAMIDINE],
    ligand_binding_affinity_index=0,
    pocket_constraints=[pocket],
    contact_constraints=[salt_bridge],
    use_potentials=True,   # required for force=True
    num_samples=3,
    folder=folder,
)
```

- `pocket_constraints`: a list of `rowan.PocketConstraint`, each placing one binder near a set of contact residues within `max_distance` angstroms.
- `contact_constraints`: a list of `rowan.ContactConstraint`, each holding two targets within `max_distance` angstroms. Setting `force=True` enforces the contact and requires `use_potentials=True`.

## Result fields

- `predictions` (also `cofolding_results`): all cofolding predictions.
- `predicted_structure_uuid` and `predicted_refined_structure_uuid`: UUIDs of the predicted structure and, when pose refinement ran, the refined one. Fetch with `rowan.retrieve_protein(uuid)`.
- `lddt`: per-residue LDDT confidence for the primary prediction.
- `scores`: confidence scores; `affinity_score`: predicted binding affinity for the primary prediction.
- `affinity_score` fields are all optional and **which ones are filled depends on the model** (the schema does not enforce this): Boltz-2 fills `pred_value`/`pred_value1`/`pred_value2` (predicted pIC50, higher = stronger; IC50 is derived from it) and `probability_binary`/`probability_binary1`/`probability_binary2` (binding probability, 0-1); Boltz-2.1 instead fills `binding_confidence` (probability the ligand is a true binder vs decoy — use for hit discovery) and `optimization_score` (binding-strength ranking, higher = stronger — use for lead optimization). The two models report non-overlapping affinity metrics, so a Boltz-2 vs Boltz-2.1 comparison shares only the `scores`.
- `posebusters_valid`: whether the primary pose passes PoseBusters.
- `strain`: ligand strain energy for the primary prediction.
