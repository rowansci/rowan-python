import rowan

# Set ROWAN_API_KEY environment variable to your API key or set rowan.api_key directly
# rowan.api_key = "rowan-sk..."

oseltamivir_SMILES = "CCOC(=O)C1=CC=CC=C1"  # ethyl benzoate for quick test

workflow = rowan.submit_admet_workflow(
    initial_smiles=oseltamivir_SMILES,
    name="Test ADMET",
)

print(f"Submitted: {workflow.uuid}")
workflow.wait_for_result().fetch_latest(in_place=True)
print(workflow)

# Test the new result() method
result = workflow.result()
print(f"{workflow=}")
print(f"{result=}")
print(f"{type(result)=}")
print(f"{result.properties=}")
