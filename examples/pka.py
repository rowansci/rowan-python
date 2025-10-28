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
    initial_molecule=Molecule.from_smiles("c1ccccc1O"),
    method="aimnet2_wagen2024",
    mode="reckless",
    name="Pyridine pKa",
)

print(result.wait_for_result().fetch_latest(in_place=True))

result2 = rowan.submit_pka_workflow(
    initial_molecule="c1ccccc1O",
    method="chemprop_nevolianis2025",
    name="Pyridine pKa (ML)",
)

print(result2.wait_for_result().fetch_latest(in_place=True))
