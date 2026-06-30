from pathlib import Path

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")
data_dir = Path(__file__).parent / "data"

protein = rowan.upload_protein("TYK2", data_dir / "tyk2_structure.pdb")
all_ligands = rowan.load_named_ligands(data_dir / "tyk2_ligands.sdf")
ligands = dict(list(all_ligands.items())[:3])

workflow = rowan.submit_binding_affinity_workflow(
    protein=protein,
    ligand_structures=list(ligands.values()),
    name="Binding Affinity — TYK2 ligands",
    folder=folder,
)
print(f"View at: https://labs.rowansci.com/binding-affinity/{workflow.uuid}")

result = workflow.result()
for name, score in zip(ligands.keys(), result.scores, strict=False):
    print(f"{name}: {score.binding_affinity:.2f} kcal/mol (strain: {score.strain})")
