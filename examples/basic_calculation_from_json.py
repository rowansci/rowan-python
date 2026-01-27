import json

import rowan

# Set ROWAN_API_KEY environment variable to your API key or set rowan.api_key directly
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
workflow.wait_for_result().fetch_latest(in_place=True)
print(workflow)
