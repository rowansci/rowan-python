from stjames import Molecule

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."

workflow = rowan.submit_basic_calculation_workflow(
    initial_molecule=Molecule.from_smiles("CC(=C)C=C"),
)

print(f"View workflow privately at: https://labs.rowansci.com/workflow/{workflow.uuid}")
result = workflow.result()
print(result)
