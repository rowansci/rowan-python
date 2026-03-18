from pathlib import Path

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."

folder = rowan.get_folder("examples")
data_dir = Path(__file__).parent / "data"

protein = rowan.upload_protein("TYK2", data_dir / "tyk2_structure.pdb")
protein.prepare()

# Step 1: build the perturbation graph from TYK2 ligands with 3D coordinates
ligands = rowan.load_named_ligands(data_dir / "tyk2_ligands.sdf")

rbfe_graph = rowan.submit_relative_binding_free_energy_graph_workflow(
    ligands=ligands,
    folder=folder,
    name="TYK2 RBFE Graph",
)
rbfe_graph_result = rbfe_graph.result()
print(rbfe_graph_result)
# <RelativeBindingFreeEnergyGraphResult edges=18 ligands=16>

# Step 2: run the FEP simulation
perturbation_workflow = rowan.submit_relative_binding_free_energy_perturbation_workflow(
    graph_result=rbfe_graph_result,
    protein=protein,
    tmd_settings="recommended",  # or "fast", "rigorous"
    folder=folder,
    name="TYK2 RBFE Perturbation",
)

print(
    f"Submitted! View at: https://labs.rowansci.com/relative-binding-free-energy-perturbation/{perturbation_workflow.uuid}"
)
print(f"UUID: {perturbation_workflow.uuid}")

# Come back later to retrieve results:
#
# workflow = rowan.retrieve_workflow(perturbation_workflow.uuid)
# result = workflow.result()
# print(result)
# # <RelativeBindingFreeEnergyPerturbationResult ligands=16>
#
# # Per-ligand ΔΔG relative to the reference
# for name, res in result.ligand_dg_results.items():
#     print(f"  {name}: ΔΔG = {res.dg:.2f} ± {res.dg_err:.2f} kcal/mol")
#
# # QC
# print(result.diagnostics)
# # RelativeBindingFreeEnergyDiagnostics(
# #     cycle_closure_rms=..., windows_completed=52, windows_failed=...
# # )
#
# # Per-edge results
# for edge in result.edges:
#     print(f"  {edge.ligand_a} → {edge.ligand_b}: ΔΔG={edge.ddg:.2f} kcal/mol")
