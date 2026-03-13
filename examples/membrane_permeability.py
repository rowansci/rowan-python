from rdkit import Chem

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
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

print("View these workflows privately:")
print(f"\thttps://labs.rowansci.com/membrane-permeability/{gnn_mtl_workflow.uuid}")
print(f"\thttps://labs.rowansci.com/mambrane-permeability/{pypermm_workflow.uuid}")
gnn_mtl_result = gnn_mtl_workflow.result()
pypermm_result = pypermm_workflow.result()

print(f"Caco-2 Papp (GNN-MTL): {gnn_mtl_result.caco2_p_app}")
print(f"Caco-2 logP (PyPermm): {pypermm_result.caco2_log_p}")
