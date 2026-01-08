from rdkit import Chem

import rowan


# rowan.api_key = ""

oseltamivir_SMILES = "C1CCOC(=O)C1=C[C@@H](OC(CC)CC)[C@H](NC(C)=O)[C@@H]([NH3+])C1CCC1"

result = rowan.submit_membrane_permeability_workflow(
    oseltamivir_SMILES,
    name="Oseltamivir Membrane Permeability",
)

print(result.wait_for_result().fetch_latest(in_place=True))

result = rowan.submit_membrane_permeability_workflow(
    Chem.MolFromSmiles(oseltamivir_SMILES),
    method="pypermm",
    name="Oseltamivir Membrane Permeability (PyPermm)",
)

print(result.wait_for_result().fetch_latest(in_place=True))
