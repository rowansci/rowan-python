"""
Run an optimization calculation on a molecule using Rowan.

See documentation at: https://docs.rowansci.com/science/quantum-chemistry/geometry-optimization
"""

from stjames import Molecule

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

workflow = rowan.submit_basic_calculation_workflow(
    initial_molecule=Molecule.from_smiles("O"),  # water
    method="GFN2-xTB",
    tasks=["optimize", "frequencies"],
    engine="xtb",
    name="Water Optimization",
    folder=folder,
)

print(f"View optimization privately at: https://labs.rowansci.com/calculation/{workflow.uuid}")
result = workflow.result()
print(result)

mols = result.molecules

print(f"\n# molecules: {len(mols)}")

print("\nAngles:")
for mol in mols:
    print(f"{mol.angle(2, 1, 3):5.2f}")

print("\nFrequencies:")
for vibrational_mode in mols[-1].vibrational_modes or []:
    print(f"{vibrational_mode.frequency:8.3f}")
