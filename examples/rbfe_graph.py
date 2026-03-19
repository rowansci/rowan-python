from pathlib import Path

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."

folder = rowan.get_folder("examples")
data_dir = Path(__file__).parent / "data"

# Load TYK2 ligands from SDF - 3D coordinates are used by the graph builder
ligands = rowan.load_named_ligands(data_dir / "tyk2_ligands.sdf")

print(f"Loaded {len(ligands)} ligands: {list(ligands.keys())}")

workflow = rowan.submit_relative_binding_free_energy_graph_workflow(
    ligands=ligands,
    folder=folder,
    name="TYK2 RBFE Graph",
)

print(
    f"View at: https://labs.rowansci.com/relative-binding-free-energy-perturbation/{workflow.uuid}"
)
result = workflow.result()
print(result)
# <RelativeBindingFreeEnergyGraphResult edges=18 ligands=16>

for edge in result.edges:
    print(f"  {edge.ligand_a} -> {edge.ligand_b}  score={edge.score:.3f}")
