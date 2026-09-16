import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

protein = rowan.create_protein_from_pdb_id("1CRN", name="crambin")
preparation_workflow = rowan.submit_protein_preparation_workflow(
    protein=protein.uuid,
    name="Prepare crambin",
    folder=folder,
)
preparation_result = preparation_workflow.result()
prepared_protein_uuid = preparation_result.prepared_protein_uuid
if prepared_protein_uuid is None:
    raise ValueError("Protein preparation returned no prepared protein")

md_workflow = rowan.submit_protein_md_workflow(
    protein=prepared_protein_uuid,
    num_trajectories=1,
    simulation_time_ns=1,
    name="MD on crambin",
    folder=folder,
)

print(f"View MD workflow privately at: https://labs.rowansci.com/protein-md/{md_workflow.uuid}")
md_result = md_workflow.result()
trajectory = md_result.trajectories[0]
print(f"Protein RMSD: {trajectory.protein_rmsd}")
print(f"Protein RMSF: {trajectory.rmsf}")
print(f"Potential energy: {trajectory.potential_energy} Hartree")
