import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

oseltamivir_SMILES = "C1CCOC(=O)C1=C[C@@H](OC(CC)CC)[C@H](NC(C)=O)[C@@H]([NH3+])C1CCC1"

workflow = rowan.submit_admet_workflow(
    initial_smiles=oseltamivir_SMILES,
    name="Oseltamivir ADMET",
    folder=folder,
)

print(f"View workflow privately at: https://labs.rowansci.com/solubility/{workflow.uuid}")
result = workflow.result()
print(result)
# e.g. <ADMETResult properties=42 preview={'mw': 180.16, 'tpsa': 75.3, ...}>
