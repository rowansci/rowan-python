"""
Calculate the redox potential of a molecule using the Rowan API.

See documentiation at: https://docs.rowansci.com/science/workflows/redox-potential
"""

from stjames import Molecule

import rowan

# Set ROWAN_API_KEY environment variable to your API key or set rowan.api_key directly
# rowan.api_key = "rowan-sk..."

workflow = rowan.submit_redox_potential_workflow(
    initial_molecule=Molecule.from_smiles("C1=CC=C(C=C1)C(=O)O"),
    mode="reckless",
    name="Benzoic Acid Redox Potential",
    oxidization=True,
    reduction=True,
)

print(f"View workflow privately at: https://labs.rowansci.com/workflow/{workflow.uuid}")
workflow.wait_for_result().fetch_latest(in_place=True)
print(workflow)
