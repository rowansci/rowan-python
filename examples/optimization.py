"""
Run an optimization calculation on a molecule using Rowan.

See documentation at: https://docs.rowansci.com/science/quantum-chemistry/geometry-optimization
"""

from stjames import Molecule

import rowan

# Set ROWAN_API_KEY environment variable to your API key or set rowan.api_key directly
# rowan.api_key = "rowan-sk..."

workflow = rowan.submit_basic_calculation_workflow(
    initial_molecule=Molecule.from_smiles("O"),
    method="GFN2-xTB",
    tasks=["optimize"],
    engine="xtb",
    name="Water Optimization",
)

print(f"View optimization privately at: https://labs.rowansci.com/workflow/{workflow.uuid}")
workflow.wait_for_result().fetch_latest(in_place=True)
print(workflow)
