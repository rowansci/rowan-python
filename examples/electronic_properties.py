"""
Calculate electronic properties (orbitals, density, ESP) using the Rowan API.

This workflow computes molecular orbitals, electron density, electrostatic potential,
atom-centered charges, bond orders, and multipole moments.

Note: This workflow does not optimize molecular geometries.

See documentation at: https://docs.rowansci.com/science/workflows/orbitals
"""

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

workflow = rowan.submit_electronic_properties_workflow(
    initial_molecule=rowan.Molecule.from_smiles("c1ccccc1"),  # benzene
    method="b97_3c",  # default: lightweight DFT
    compute_density_cube=True,
    compute_electrostatic_potential_cube=True,
    compute_num_occupied_orbitals=3,  # HOMO, HOMO-1, HOMO-2
    compute_num_virtual_orbitals=3,  # LUMO, LUMO+1, LUMO+2
    name="Benzene electronic properties",
    folder=folder,
)

print(f"View workflow privately at: https://labs.rowansci.com/orbitals/{workflow.uuid}")
result = workflow.result()
print(result)
# e.g. <ElectronicPropertiesResult dipole=1.85 D>
