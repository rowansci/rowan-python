import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."

oseltamivir_SMILES = "C1CCOC(=O)C1=C[C@@H](OC(CC)CC)[C@H](NC(C)=O)[C@@H]([NH3+])C1CCC1"

workflow = rowan.submit_solubility_workflow(
    initial_smiles=oseltamivir_SMILES,
    solubility_method="kingfisher",
    solvents=["O"],
    temperatures=[298.15],
    name="Oseltamivir aqueous solubility",
)

print(f"View workflow privately at: https://labs.rowansci.com/workflow/{workflow.uuid}")
result = workflow.result()
print(result)  # Solubility in log(mol/L), temperature in kelvin
