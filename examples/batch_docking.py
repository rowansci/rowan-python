import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

ligands = [
    "CCC(C)(C)NC1=NCC2(CCC(=O)C2C)N1",
    "CCC(C)CN=C1NCC2(CCCOC2)CN1",
    "CC(C)CCNC1=NCC2CC(COC2=N)O1",
    "CCC(CC)NC1=NCC2CC(CO)CC12",
    "CCC(C)CN=C1NC=C2CCC(O)CC2=N1",
]

protein = rowan.create_protein_from_pdb_id(
    "CDK2", "1HCK", project_uuid=rowan.default_project().uuid
)

protein.sanitize()

workflow = rowan.submit_batch_docking_workflow(
    ligands,
    protein.uuid,
    pocket=[[103.55, 100.59, 82.99], [27.76, 32.67, 48.79]],
    executable="qvina2",
    scoring_function="vina",
    folder=folder,
)


print(f"View workflow privately at: https://labs.rowansci.com/batch-docking/{workflow.uuid}")
print(f"Workflow UUID: {workflow.uuid}")

# Stream partial scores as each ligand completes; final iteration is the complete result.
for result in workflow.stream_result(poll_interval=30):
    completed = sum(s is not None for s in result.scores.values())
    print(f"  {completed}/{len(ligands)} complete...")

print(result.scores)  # dict of SMILES → best docking score
