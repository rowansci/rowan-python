"""
Calculate the redox potential of a molecule using the Rowan API.

See documentiation at: https://docs.rowansci.com/science/workflows/redox-potential
"""

from stjames import Molecule

import rowan

# rowan.api_key = ""

result = rowan.submit_redox_potential_workflow(
    initial_molecule=Molecule.from_smiles("C1=CC=C(C=C1)C(=O)O"),
    mode="reckless",
    name="Benzoic Acid Redox Potential",
    oxidization=True,
    reduction=True,
)

print(result.wait_for_result().fetch_latest(in_place=True))
