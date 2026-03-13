import json

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."

with open("examples/data/workflow_example.json") as f:
    workflow_data = json.load(f)

workflow = rowan.submit_workflow(
    workflow_type="basic_calculation",
    workflow_data=workflow_data,
    name="basic calculation from json",
    initial_molecule=workflow_data["initial_molecule"],
)

print(f"View workflow privately at: https://labs.rowansci.com/calculation/{workflow.uuid}")
result = workflow.result()
print(result)
