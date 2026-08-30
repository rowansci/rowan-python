import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

# Dasatinib — redocked into its own ABL1 co-crystal structure (PDB: 2GQG)
dasatinib = rowan.Molecule.from_smiles("Cc1nc(Nc2ncc(C(=O)Nc3c(C)cccc3Cl)s2)cc(N2CCN(CCO)CC2)n1")

protein = rowan.create_protein_from_pdb_id("2GQG")
protein = protein.select_chains(["A"])
preparation_workflow = rowan.submit_protein_preparation_workflow(
    protein=protein.uuid,
    add_missing_method="pdbfixer",
    name="Prepare ABL1",
    folder=folder,
)
prepared_protein_uuid = preparation_workflow.result().prepared_protein_uuid

# Pocket is [[center_x, center_y, center_z], [size_x, size_y, size_z]] in Å.
center = [44.59, 79.75, 39.59]
size = [24.15, 21.33, 19.88]

# Induced-fit docking soft-docks candidate poses, relaxes the receptor around each with
# restrained local minimization, and redocks into the relaxed receptor. It requires
# VinaSettings with executable="vina" or "qvina2" (the default is "vina").
induced_fit_settings = rowan.InducedFitSettings(
    max_receptors=6,
    flexible_sidechain_radius=5.0,
)

workflow = rowan.submit_docking_workflow(
    prepared_protein_uuid,
    pocket=[center, size],
    initial_molecule=dasatinib,
    induced_fit_settings=induced_fit_settings,
    name="Dasatinib induced-fit docking",
    folder=folder,
)

print(f"View workflow privately at: https://labs.rowansci.com/docking/{workflow.uuid}")

result = workflow.result()
print(result)

# With induced-fit docking, poses are ranked by induced_fit_score rather than raw docking score.
for i, score in enumerate(result.scores):
    print(
        f"  Pose {i}: induced_fit_score={score.induced_fit_score:.3f}  "
        f"score={score.score:.3f}  receptor_strain={score.receptor_strain}"
    )

# Download the top-scoring complex and best-ranked relaxed receptor as PDBs
complex_protein = result.get_complex(0)
complex_protein.download_pdb_file(name="dasatinib_2GQG_induced_complex")
print("Saved dasatinib_2GQG_induced_complex.pdb")

induced_receptors = result.get_induced_receptors()
if induced_receptors:
    induced_receptors[0].download_pdb_file(name="dasatinib_2GQG_induced_receptor")
    print("Saved dasatinib_2GQG_induced_receptor.pdb")
