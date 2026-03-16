import stjames

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples/docking-screen")

HARTREE_TO_KCALMOL = 627.5096

ligands = [
    "CCC(C)(C)NC1=NCC2(CCC(=O)C2C)N1",
    "CCC(C)CN=C1NCC2(CCCOC2)CN1",
    "CC(C)CCNC1=NCC2CC(COC2=N)O1",
    "CCC(CC)NC1=NCC2CC(CO)CC12",
    "CCC(C)CN=C1NC=C2CCC(O)CC2=N1",
]


workflows = []
results = {}

protein = rowan.create_protein_from_pdb_id(
    "CDK2", "1HCK", project_uuid=rowan.default_project().uuid
)

protein.sanitize()

for ligand in ligands:
    workflow = rowan.submit_docking_workflow(
        protein.uuid,
        pocket=[[103.55, 100.59, 82.99], [27.76, 32.67, 48.79]],
        initial_molecule=stjames.Molecule.from_smiles(ligand),
        name=f"Docking {ligand}",
        folder=folder,
    )
    workflows.append(workflow)


workflow_results = [(w, w.result()) for w in workflows]

lowest_conformer_energy = 0
for workflow, result in workflow_results:
    for conformer_uuid in result.conformers:
        energy = rowan.retrieve_calculation_molecules(conformer_uuid)[0]["energy"]
        lowest_conformer_energy = min(lowest_conformer_energy, energy)

    sorted_scores = sorted(result.scores, key=lambda s: s.score)
    for score in sorted_scores:
        pose_energy = rowan.retrieve_calculation_molecules(score.pose)[0]["energy"]
        strain = (pose_energy - lowest_conformer_energy) * HARTREE_TO_KCALMOL
        if score.posebusters_valid and strain < 4:
            results[workflow.name] = score
            break

print(results)
