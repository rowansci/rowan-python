"""Detect protein binding pockets with mixed-solvent molecular dynamics."""

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
if prepared_protein_uuid is None:
    raise ValueError("Protein preparation returned no prepared protein")

# Runs three independent replicas, each with 5 ns equilibration and 15 ns production.
workflow = rowan.submit_pocket_detection_workflow(
    protein=prepared_protein_uuid,
    settings=rowan.MixMDSettings(probe=rowan.Solvent.ACETONITRILE),
    name="MixMD pocket detection on thymidine phosphorylase",
    folder=folder,
)
print(f"View workflow privately at: https://labs.rowansci.com/pocket-detection/{workflow.uuid}")
result = workflow.result()

for message in result.messages:
    print(message)

print(f"Detected {len(result.pockets)} pocket(s):")
for i, pocket in enumerate(result.pockets):
    print(f"  Pocket {i}: score={pocket.score}, enriched-region volume={pocket.volume} Å³")
    print(f"    center={pocket.pocket_center}, sides={pocket.pocket_sides}")
    print(f"    residues={pocket.residue_numbers}")

# Pocket coordinates refer to this reference structure.
reference_protein = result.get_reference_protein()
if reference_protein is not None:
    print(f"Reference protein: {reference_protein.uuid}")
