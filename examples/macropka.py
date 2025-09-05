import rowan

# Set ROWAN_API_KEY environment variable to your API key or set rowan.api_key directly
# rowan.api_key = "rowan-sk..."

result = rowan.submit_macropka_workflow(
    # oseltamivir
    initial_smiles="C1CCOC(=O)C1=C[C@@H](OC(CC)CC)[C@H](NC(C)=O)[C@@H]([NH3+])C1CCC1",
    name="Oseltamivir macropka",
    compute_aqueous_solubility=True,
)


# Solubility units are log(mol/L)
print(result.wait_for_result().fetch_latest(in_place=True))
