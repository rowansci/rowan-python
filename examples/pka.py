"""
Calculate the pKa of the conjugate acid of pyridine using the Rowan API.

Experimental value ≈5.23

See documentiation at: https://docs.rowansci.com/science/workflows/pka
and preprint at: https://rowansci.com/publications/pka-prediction
"""

from stjames import Molecule

import rowan

# rowan.api_key = ""

result = rowan.submit_pka_workflow(
    initial_molecule=Molecule.from_smiles("n1ccccc1"),
    mode="reckless",
    name="Pyridine pKa",
)

print(result.wait_for_result().fetch_latest(in_place=True))
