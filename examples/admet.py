import rowan

# rowan.api_key = ""

oseltamivir_SMILES = "C1CCOC(=O)C1=C[C@@H](OC(CC)CC)[C@H](NC(C)=O)[C@@H]([NH3+])C1CCC1"

result = rowan.submit_admet_workflow(
    initial_smiles=oseltamivir_SMILES,
    name="Oseltamivir ADMET",
)

print(result.wait_for_result().fetch_latest(in_place=True))
