import stjames
from stjames import Molecule

import rowan

# Set ROWAN_API_KEY environment variable to your API key or set rowan.api_key directly
# rowan.api_key = "rowan-sk..."

conformer_dependent_redox_folder = rowan.create_folder(name="Conformer dependent redox workflows")

workflow = rowan.submit_conformer_search_workflow(
    initial_molecule=Molecule.from_smiles("CC(C)Cc1ccc(C(=O)c2ccc(O)cc2)cc1"),
    folder_uuid=conformer_dependent_redox_folder.uuid,
)
print(f"View workflow privately at: https://labs.rowansci.com/workflow/{workflow.uuid}")
workflow.wait_for_result().fetch_latest(in_place=True)

redox_potential_workflows = []

for conformer in workflow.data["conformer_uuids"][:10]:
    uuid = conformer[0]
    molecule = rowan.retrieve_calculation_molecules(uuid)[0]
    stjames_molecule = stjames.Molecule.model_validate(molecule)
    redox_potential_workflows.append(
        rowan.submit_redox_potential_workflow(
            stjames_molecule,
            reduction=True,
            oxidization=True,
            folder_uuid=conformer_dependent_redox_folder.uuid,
        )
    )

for workflow in redox_potential_workflows:
    workflow.wait_for_workflow()
    workflow.fetch_latest(in_place=True)


print([workflow.data["redox_potential"] for workflow in redox_potential_workflows])
