"""
Retrieve a previously run water optimization and frequency calculation.
"""

from stjames import Molecule

import rowan

# Set ROWAN_API_KEY environment variable to your API key or set rowan.api_key directly
# rowan.api_key = "rowan-sk..."

workflow_uuid = "2e8632bf-6d68-4099-b0c7-96378c7b7261"
workflow = rowan.retrieve_workflow(workflow_uuid)
print(f"View workflow at: https://labs.rowansci.com/calculation/{workflow_uuid}")

calc_uuid = workflow.data["calculation_uuid"]
mol_json = rowan.retrieve_calculation_molecules(calc_uuid, return_frequencies=True)
mols = [Molecule.model_validate(json) for json in mol_json]

print(f"\n# molecules: {len(mols)}")

print("\nAngles:")
for mol in mols:
    print(f"{mol.angle(2, 1, 3):5.2f}")

print("\nFrequencies:")
for vibrational_mode in mols[-1].vibrational_modes:
    print(f"{vibrational_mode.frequency:8.3f}")
