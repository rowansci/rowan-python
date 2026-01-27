from stjames import Molecule

import rowan

# Set ROWAN_API_KEY environment variable to your API key or set rowan.api_key directly
# rowan.api_key = "rowan-sk..."

workflow = rowan.submit_basic_calculation_workflow(
    initial_molecule=Molecule.from_smiles("CC(=C)C=C"),
)

print(f"View workflow privately at: https://labs.rowansci.com/workflow/{workflow.uuid}")
workflow.wait_for_result().fetch_latest(in_place=True)
print(workflow)
