import time

import rowan

# Set ROWAN_API_KEY environment variable to your API key or set rowan.api_key directly
# rowan.api_key = "rowan-sk..."

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
time.sleep(60)
protein.refresh()

workflow = rowan.submit_batch_docking_workflow(
    ligands,
    protein.uuid,
    pocket=[[103.55, 100.59, 82.99], [27.76, 32.67, 48.79]],
    executable="qvina2",
    scoring_function="vina",
)


print(f"View workflow privately at: https://labs.rowansci.com/batch-docking/{workflow.uuid}")
workflow.wait_for_result().fetch_latest(in_place=True)
print(workflow.data["best_scores"])
