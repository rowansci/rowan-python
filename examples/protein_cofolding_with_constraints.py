"""Cofolding with both pocket and contact constraints.

Example: bovine trypsin + benzamidine, a textbook S1-pocket binder.
Trypsin's S1 specificity pocket is anchored by Asp189 at its base, and
benzamidine binds via a bidentate salt bridge between its amidinium group
and the Asp189 carboxylate.

We encode this prior knowledge with two constraints:

1. PocketConstraint - the ligand must sit inside the S1 pocket, near the
   catalytic triad (His57, Asp102, Ser195) and Asp189.
2. ContactConstraint - the amidine carbon of benzamidine must come within
   ~4 A of Asp189 (the salt bridge interaction).

The four residues are located by their conserved chymotrypsin-family motifs
in the sequence below: AAHCY (His), NNDIML (Asp triad), SDSSCK (S1 Asp
before Cys191), GDSGGP (Ser).

ConstraintTargets are constructed explicitly with (input_type, input_index,
token_index) so the addressing scheme is unambiguous: input_index selects which
protein / ligand in the input lists, and token_index is the residue index for
proteins and the atom index for ligands. The same three fields appear in the
cofolding submit form on labs.rowansci.com.
"""

import rowan

# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

TRYPSIN = (
    "IVGGYTCGANTVPYQVSLNSGYHFCGGSLINSQWVVSAAHCYKSGIQVRLGEDNINVVEGNEQFISASKSIVHPSYNSNTLNN"
    "DIMLIKLKSAAYTSYDVPLGTQCLISGWGNTKSSGTSYPDVLKCLKAPILSDSSCKSAYPGQITSNMFCAGYLEGGKDSCQGD"
    "SGGPVVCSGKLQGIVSWGSGCAQKNKPGVYTKVCNYVSWIKQTIASN"
)
BENZAMIDINE = "NC(=N)c1ccccc1"

# Targets for the four pocket-defining residues on the trypsin chain (protein 0).
his57 = rowan.ConstraintTarget(input_type="protein", input_index=0, token_index=39)  # AAHCY
asp102 = rowan.ConstraintTarget(input_type="protein", input_index=0, token_index=83)  # NNDIML
asp189 = rowan.ConstraintTarget(input_type="protein", input_index=0, token_index=134)  # SDSSCK
ser195 = rowan.ConstraintTarget(input_type="protein", input_index=0, token_index=166)  # GDSGGP

# Target for the central amidine carbon of benzamidine (atom 1 in NC(=N)c1ccccc1).
benzamidine_amidine_c = rowan.ConstraintTarget(input_type="ligand", input_index=0, token_index=1)

pocket = rowan.PocketConstraint(
    input_type="ligand",  # the binder being placed into the pocket
    input_index=0,  # benzamidine is the first (and only) ligand
    contacts=[his57, asp102, asp189, ser195],
    max_distance=6.0,
)

salt_bridge = rowan.ContactConstraint(
    token_1=benzamidine_amidine_c,
    token_2=asp189,
    max_distance=4.0,
    force=True,
)

workflow = rowan.submit_protein_cofolding_workflow(
    initial_protein_sequences=[TRYPSIN],
    initial_smiles_list=[BENZAMIDINE],
    ligand_binding_affinity_index=0,
    pocket_constraints=[pocket],
    contact_constraints=[salt_bridge],
    use_potentials=True,  # required for force=True
    num_samples=3,
    name="Trypsin + Benzamidine (constrained)",
    do_pose_refinement=True,
    folder=folder,
)

print(f"View workflow privately at: https://labs.rowansci.com/protein-cofolding/{workflow.uuid}")
result = workflow.result()
print(result)
for i, pred in enumerate(result.predictions):
    iptm = pred.scores.iptm if pred.scores else None
    print(f"  sample {i}: iptm={iptm}")
