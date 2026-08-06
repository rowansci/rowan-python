import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

protein = rowan.create_protein_from_pdb_id("1CRN", name="crambin")

workflow = rowan.submit_protein_preparation_workflow(
    protein=protein.uuid,
    name="Prepare crambin",
    folder=folder,
)

print(f"View workflow privately at: https://labs.rowansci.com/protein-preparation/{workflow.uuid}")
result = workflow.result()
print(result)
print(result.prepared_protein_uuid)
