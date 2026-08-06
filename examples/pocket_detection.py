import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

protein = rowan.create_protein_from_pdb_id("1OTP", name="thymidine phosphorylase")
preparation_workflow = rowan.submit_protein_preparation_workflow(
    protein=protein.uuid,
    name="Prepare thymidine phosphorylase",
    folder=folder,
)
prepared_protein_uuid = preparation_workflow.result().prepared_protein_uuid

workflow = rowan.submit_pocket_detection_workflow(
    protein=prepared_protein_uuid,
    name="Pocket detection on thymidine phosphorylase",
    folder=folder,
)

result = workflow.result()

print(f"Detected {len(result.pockets)} pocket(s):")
for i, pocket in enumerate(result.pockets):
    print(f"  Pocket {i}: score={pocket.score}, volume={pocket.volume} Å³")
    print(f"    center={pocket.pocket_center}")
    print(f"    sides={pocket.pocket_sides}")
    print(f"    residues={pocket.residue_numbers}")

print(f"View workflow privately at: https://labs.rowansci.com/pocket-detection/{workflow.uuid}")
