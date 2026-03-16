import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

protein = rowan.create_protein_from_pdb_id(
    "crambin", "1CRN", project_uuid=rowan.default_project().uuid
)

protein.sanitize()


md_workflow = rowan.submit_protein_md_workflow(
    protein=protein,
    num_trajectories=1,
    simulation_time_ns=1,
    name="MD on crambin",
    folder=folder,
)

print(f"View MD workflow privately at: https://labs.rowansci.com/protein-md/{md_workflow.uuid}")
