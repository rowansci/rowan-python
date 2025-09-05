import rowan

# rowan.api_key = ""


oseltamivir_SMILES = "C1CCOC(=O)C1=C[C@@H](OC(CC)CC)[C@H](NC(C)=O)[C@@H]([NH3+])C1CCC1"

result = rowan.submit_solubility_workflow(
    initial_smiles=oseltamivir_SMILES,
    solubility_method="kingfisher",
    solvents=["O"],
    temperatures=[298.15],
    name="Oseltamivir aqueous solubility",
)

# Solubility units are log(mol/L) and temperature units are kelvin
print(result.wait_for_result().fetch_latest(in_place=True))
