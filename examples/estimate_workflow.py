"""Compare dispatch info for different property predictions to decide which to submit."""

import rowan

folder = rowan.get_folder("examples")
smiles = "c1ccc(CC(=O)O)cc1"  # phenylacetic acid

# Create drafts, nothing runs yet
admet_draft = rowan.submit_admet_workflow(
    initial_smiles=smiles, name="Estimate Workflow - ADMET", folder=folder, is_draft=True
)
solubility_draft = rowan.submit_solubility_workflow(
    initial_smiles=smiles, name="Estimate Workflow - Solubility", folder=folder, is_draft=True
)
membrane_draft = rowan.submit_membrane_permeability_workflow(
    initial_molecule=smiles,
    name="Estimate Workflow - Membrane Permeability",
    folder=folder,
    is_draft=True,
)

# Check dispatch info for each
print("ADMET:", admet_draft.dispatch_info())
print("Solubility:", solubility_draft.dispatch_info())
print("Membrane:", membrane_draft.dispatch_info())

# Pick ADMET, clean up the rest
solubility_draft.delete()
membrane_draft.delete()
admet_draft.submit_draft()

result = admet_draft.result()
print(result)
