import rowan

folder = rowan.get_folder("examples")

tg2_inhibitor = rowan.Molecule.from_smiles(
    "CCC(=O)NCCCC[C@H](NC(=O)Cc1ccc(Cl)c(Cl)c1)C(=O)N1CCN(C(=O)c2cccc3ccccc23)CC1"
)

protein = rowan.create_protein_from_pdb_id("2Q3Z")
protein = protein.select_chains(["A"])
preparation_workflow = rowan.submit_protein_preparation_workflow(
    protein=protein.uuid,
    name="Prepare TG2",
    folder=folder,
)
protein = preparation_workflow.result().get_prepared_protein()

center = [-1.079, -3.081, 18.122]
size = [22.22, 14.08, 21.74]

# Protein preparation renumbers Cys277 to residue 278.
cys277_sg_index = protein.get_atom_index(chain="A", residue=278, atom="SG")
gnina_settings = rowan.GninaSettings(
    scoring_function="vina",
    covalent_ligand_atom_index=0,
    covalent_protein_atom_index=cys277_sg_index,
)

workflow = rowan.submit_docking_workflow(
    protein.uuid,
    pocket=[center, size],
    initial_molecule=tg2_inhibitor,
    docking_settings=gnina_settings,
    name="TG2 covalent docking (Cys277, 2Q3Z)",
    folder=folder,
)

print(f"View workflow privately at: https://labs.rowansci.com/docking/{workflow.uuid}")

result = workflow.result()
print(result)

for i, score in enumerate(result.scores):
    print(f"  Pose {i}: score={score.score:.3f}")
