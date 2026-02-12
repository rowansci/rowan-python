"""
Perform a multistage geometry optimization using the Rowan API.

Available modes:
- "reckless": Fastest, least accurate
- "rapid": Good balance of speed and accuracy (default)
- "careful": More accurate, slower
- "meticulous": Most accurate, slowest

See documentation at: https://docs.rowansci.com/science/workflows/multistage-optimization
"""

from stjames import Molecule

import rowan

# Set ROWAN_API_KEY environment variable to your API key or set rowan.api_key directly
# rowan.api_key = "rowan-sk..."

workflow = rowan.submit_multistage_optimization_workflow(
    initial_molecule=Molecule.from_smiles("C1CCC1"),  # cyclobutane
    mode="rapid",
    name="Multistage optimization cyclobutane",
)

print(f"View workflow privately at: https://labs.rowansci.com/multistage-opt/{workflow.uuid}")
workflow.wait_for_result().fetch_latest(in_place=True)
print(workflow.data)
