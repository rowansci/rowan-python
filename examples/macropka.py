import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

workflow = rowan.submit_macropka_workflow(
    # oseltamivir
    initial_smiles="C1CCOC(=O)C1=C[C@@H](OC(CC)CC)[C@H](NC(C)=O)[C@@H]([NH3+])C1CCC1",
    name="Oseltamivir macropka",
    compute_aqueous_solubility=True,
    folder=folder,
)

print(f"View workflow privately at: https://labs.rowansci.com/macropka/{workflow.uuid}")
result = workflow.result()
print(result)
# e.g. <MacropKaResult isoelectric_point=6.5 pka_values=3>
