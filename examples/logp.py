import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

oseltamivir_SMILES = "C1CCOC(=O)C1=C[C@@H](OC(CC)CC)[C@H](NC(C)=O)[C@@H]([NH3+])C1CCC1"

workflow = rowan.submit_logp_workflow(
    initial_smiles=oseltamivir_SMILES,
    method="chemprop_sangster2026",
    name="Oseltamivir logP",
    folder=folder,
)

print(f"View workflow privately at: https://labs.rowansci.com/logp/{workflow.uuid}")
result = workflow.result()
print(result)
# e.g. <LogPResult logp=2.541>
