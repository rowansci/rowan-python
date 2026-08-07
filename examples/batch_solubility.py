import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

smileses = ["CCO", "CC(=O)O", "c1ccccc1"]
workflows = rowan.submit_solubility_workflow_group(
    initial_smileses=smileses,
    method="kingfisher",
    solvents=["water"],
    temperatures=[298.15],
    names=["Ethanol", "Acetic acid", "Benzene"],
    folder=folder,
)

for workflow in workflows:
    print(f"View workflow privately at: https://labs.rowansci.com/solubility/{workflow.uuid}")

results = [workflow.result() for workflow in workflows]
print(results)
