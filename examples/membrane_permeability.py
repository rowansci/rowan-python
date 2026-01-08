import rowan

from rdkit import Chem

smiles = "CC1=C(N=CN1)CSCCNC(=NC)NC#N"

gnn_mtl_result = rowan.submit_membrane_permeability_workflow(
    smiles,
    name="GNN=MTL permeability",
)

pypermm_result = rowan.submit_membrane_permeability_workflow(
    Chem.MolFromSmiles(smiles),
    method="pypermm",
    name="Oseltamivir Membrane Permeability (PyPermm)",
)

gnn_mtl_result.wait_for_result().fetch_latest(in_place=True)
pypermm_result.wait_for_result().fetch_latest(in_place=True)


print(f"Caco-2 Papp (GNN-MTL): {gnn_mtl_result.data['caco_2_P_app']}")
print(f"Caco-2 P0 (PyPermm): {pypermm_result.data['caco_2_logP']}")
