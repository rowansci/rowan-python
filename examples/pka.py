"""
Calculate the pKa of the conjugate acid of pyridine using the Rowan API.

Experimental value ≈5.23

See documentiation at: https://docs.rowansci.com/science/workflows/pka
and preprint at: https://rowansci.com/publications/pka-prediction
"""

from stjames import Molecule

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples/pka")

workflow = rowan.submit_pka_workflow(
    initial_molecule=Molecule.from_smiles("c1ccccc1O"),
    method="aimnet2_wagen2024",
    mode="reckless",
    name="Pyridine pKa",
    folder=folder,
)

print(
    f"View pKa with aimnet2_wagen2024 privately at: https://labs.rowansci.com/pka/{workflow.uuid}"
)
result = workflow.result()
print(result)
# e.g. <pKaResult acid=4.75>

workflow2 = rowan.submit_pka_workflow(
    initial_molecule="c1ccccc1O",
    method="chemprop_nevolianis2025",
    name="Pyridine pKa (ML)",
    folder=folder,
)

print(
    f"View pKa with chemprop_nevolianis2025 privately at: https://labs.rowansci.com/pka/{workflow2.uuid}"
)
result2 = workflow2.result()
print(result2)
# e.g. <pKaResult acid=4.75>
