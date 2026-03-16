"""
Perform a multistage geometry optimization using the Rowan API.

Available modes:
- "reckless": Fastest, least accurate
- "rapid": Good balance of speed and accuracy (default)
- "careful": More accurate, slower
- "meticulous": Most accurate, slowest

See documentation at: https://docs.rowansci.com/science/workflows/multistage-optimization
"""

import time

from stjames import Molecule

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

workflow = rowan.submit_multistage_optimization_workflow(
    initial_molecule=Molecule.from_smiles("C1CCC1"),  # cyclobutane
    mode="rapid",
    name="Multistage optimization cyclobutane",
    folder=folder,
)

print(f"View workflow privately at: https://labs.rowansci.com/multistage-opt/{workflow.uuid}")
print(f"Workflow UUID: {workflow.uuid}")

# Poll stage-by-stage progress while running.
while not workflow.done():
    partial = workflow.result(wait=False)
    print(f"  {len(partial.calculation_uuids)} stages done, energy={partial.energy}")
    time.sleep(10)

result = workflow.result()
print(result)

# Or pick up a completed workflow later by UUID:
# result = rowan.retrieve_workflow(workflow.uuid).result(wait=False)
