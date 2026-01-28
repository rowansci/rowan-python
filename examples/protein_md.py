import time

import rowan

protein = rowan.create_protein_from_pdb_id(
    "crambin", "1CRN", project_uuid=rowan.default_project().uuid
)

protein.sanitize()
time.sleep(60)
protein.refresh()


md_workflow = rowan.submit_protein_md_workflow(
    protein=protein,
    num_trajectories=1,
    simulation_time_ns=1,
    name="MD on crambin",
)

print(f"View MD workflow privately at: https://labs.rowansci.com/protein-md/{md_workflow.uuid}")
