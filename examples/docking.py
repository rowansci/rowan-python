import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

# Dasatinib — redocked into its own ABL1 co-crystal structure (PDB: 2GQG)
dasatinib = rowan.Molecule.from_smiles("Cc1nc(Nc2ncc(C(=O)Nc3c(C)cccc3Cl)s2)cc(N2CCN(CCO)CC2)n1")

protein = rowan.create_protein_from_pdb_id("2GQG")
if len(protein.chains) > 1:
    protein = protein.select_chains([protein.chains[0]])
protein.prepare()

# Pocket is [[center_x, center_y, center_z], [size_x, size_y, size_z]] in Å.
# For a co-crystal structure, extract these from the bound ligand's position.
# Use rowan.submit_pocket_detection_workflow() when the binding site is unknown.
center = [44.59, 79.75, 39.59]
size = [24.15, 21.33, 19.88]
workflow = rowan.submit_docking_workflow(
    protein,
    pocket=[center, size],
    initial_molecule=dasatinib,
    name="dasatinib docking (2GQG redock)",
    folder=folder,
)

print(f"View workflow privately at: https://labs.rowansci.com/docking/{workflow.uuid}")

result = workflow.result()
print(result)

for i, score in enumerate(result.scores):
    print(f"  Pose {i}: score={score.score:.3f}  posebusters_valid={score.posebusters_valid}")

# Download the top-scoring protein–ligand complex as a PDB
complex_protein = result.get_complex(0)
complex_protein.download_pdb_file("dasatinib_2GQG_complex.pdb")
print("Saved dasatinib_2GQG_complex.pdb")
