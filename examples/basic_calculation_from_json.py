import json

import rowan

# rowan.api_key = ""

with open("examples/data/workflow_example.json") as f:
    workflow_data = json.load(f)

result = (
    rowan.submit_workflow(
        workflow_type="basic_calculation",
        workflow_data=workflow_data,
        name="basic calculation from json",
        initial_molecule=workflow_data["initial_molecule"],
    )
    .wait_for_result()
    .fetch_latest()
)

print(result)
