import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

workflow = rowan.submit_tautomer_search_workflow(
    initial_molecule=rowan.Molecule.from_smiles("O=C1C=CC=CN1"),
    mode="reckless",
    name="2-Pyridone Tautomers",
    folder=folder,
)

print(f"View workflow privately at: https://labs.rowansci.com/tautomers/{workflow.uuid}")
result = workflow.result()
print(result)
# e.g. <TautomerResult tautomers=3 lowest_energy=-323.456789 H>
