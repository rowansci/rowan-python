"""
Calculate the pKa of the conjugate acid of pyridine using the Rowan API.

Experimental value ≈5.23

See documentiation at: https://docs.rowansci.com/science/workflows/pka
and preprint at: https://rowansci.com/publications/pka-prediction
"""

from stjames import Molecule

import rowan

# Set ROWAN_API_KEY environment variable to your API key or set rowan.api_key directly
# rowan.api_key = "rowan-sk..."

workflow = rowan.submit_pka_workflow(
    initial_molecule=Molecule.from_smiles("c1ccccc1O"),
    method="aimnet2_wagen2024",
    mode="reckless",
    name="Pyridine pKa",
)

print(
    f"View pKa with aimnet2_wagen2024 privately at: https://labs.rowansci.com/pka/{workflow.uuid}"
)
workflow.wait_for_result().fetch_latest(in_place=True)
print(workflow)

workflow2 = rowan.submit_pka_workflow(
    initial_molecule="c1ccccc1O",
    method="chemprop_nevolianis2025",
    name="Pyridine pKa (ML)",
)

print(
    f"View pKa with chemprop_nevolianis2025 privately at: https://labs.rowansci.com/pka/{workflow2.uuid}"
)
workflow2.wait_for_result().fetch_latest(in_place=True)
print(workflow2)
