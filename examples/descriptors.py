"""
Calculate molecular descriptors using the Rowan API.

Computes molecular descriptors including COSMO descriptors (surface area,
screening charge, dielectric energy, polar surface area) in water by default.

See documentation at: https://docs.rowansci.com/science/workflows/descriptors
"""

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

workflow = rowan.submit_descriptors_workflow(
    initial_molecule=rowan.Molecule.from_smiles("CC(=O)Oc1ccccc1C(=O)O"),
    name="Aspirin Descriptors",
    folder=folder,
)

print(f"View workflow privately at: https://labs.rowansci.com/descriptors/{workflow.uuid}")
result = workflow.result()
print(result)
# e.g. <DescriptorsResult n=1686>
