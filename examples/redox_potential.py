"""
Calculate the redox potential of a molecule using the Rowan API.

See documentiation at: https://docs.rowansci.com/science/workflows/redox-potential
"""

from stjames import Molecule

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."

workflow = rowan.submit_redox_potential_workflow(
    initial_molecule=Molecule.from_smiles("C1=CC=C(C=C1)C(=O)O"),
    mode="reckless",
    name="Benzoic Acid Redox Potential",
    oxidization=True,
    reduction=True,
)

print(f"View workflow privately at: https://labs.rowansci.com/redox-potential/{workflow.uuid}")
result = workflow.result()
print(result)
