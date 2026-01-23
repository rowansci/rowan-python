from rdkit import Chem

import rowan

# Set ROWAN_API_KEY environment variable to your API key or set rowan.api_key directly
# rowan.api_key = "rowan-sk..."

smiles = "CC1=C(N=CN1)CSCCNC(=NC)NC#N"

gnn_mtl_workflow = rowan.submit_membrane_permeability_workflow(
    smiles,
    name="GNN=MTL permeability",
)

pypermm_workflow = rowan.submit_membrane_permeability_workflow(
    Chem.MolFromSmiles(smiles),
    method="pypermm",
    name="Oseltamivir Membrane Permeability (PyPermm)",
)

print(
    f"View GNN=MTL permeability privately at: https://labs.rowansci.com/workflow/{gnn_mtl_workflow.uuid}"
)
print(
    f"View Oseltamivir membrane permeability privately at: https://labs.rowansci.com/workflow/{pypermm_workflow.uuid}"
)
gnn_mtl_workflow.wait_for_result().fetch_latest(in_place=True)
pypermm_workflow.wait_for_result().fetch_latest(in_place=True)

print(f"Caco-2 Papp (GNN-MTL): {gnn_mtl_workflow.data['caco_2_P_app']}")
print(f"Caco-2 P0 (PyPermm): {pypermm_workflow.data['caco_2_logP']}")
print("View these workflows privately:")
print(f"\thttps://labs.rowansci.com/workflow/{gnn_mtl_workflow.uuid}")
print(f"\thttps://labs.rowansci.com/workflow/{pypermm_workflow.uuid}")
