"""
Calculate the redox potential of a molecule using the Rowan API.

See documentiation at: https://docs.rowansci.com/science/workflows/redox-potential
"""

from stjames import Molecule

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

workflow = rowan.submit_redox_potential_workflow(
    initial_molecule=Molecule.from_smiles("C1=CC=C(C=C1)C(=O)O"),
    mode="reckless",
    name="Benzoic Acid Redox Potential",
    oxidation=True,
    reduction=True,
    folder=folder,
)

print(f"View workflow privately at: https://labs.rowansci.com/redox-potential/{workflow.uuid}")
result = workflow.result()
print(result)
# e.g. <RedoxPotentialResult oxidation=1.234V reduction=-0.567V>
