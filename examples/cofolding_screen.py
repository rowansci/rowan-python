import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples/cofolding-screen")

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

for ligand in ligands:
    workflow = rowan.submit_protein_cofolding_workflow(
        initial_protein_sequences=[
            "MENFQKVEKIGEGTYGVVYKARNKLTGEVVALKKIRLDTETEGVPSTAIREISLLKELNHPNIVKLLDVIHTENKLYLVFEFLHQDLKKFMDASALTGIPLPLIKSYLFQLLQGLAFCHSHRVLHRDLKPQNLLINTEGAIKLADFGLARAFGVPVRTYTHEVVTLWYRAPEILLGCKYYSTAVDIWSLGCIFAEMVTRRALFPGDSEIDQLFRIFRTLGTPDEVVWPGVTSMPDYKPSFPKWARQDFSKVVPPLDEDGRSLLSQMLHYDPNKRISAKAALAHPFFQDVTKPVPHLRL"
        ],
        initial_smiles_list=[ligand],
        ligand_binding_affinity_index=0,
        name=f"Cofolding {ligand}",
        folder=folder,
    )
    workflows.append(workflow)

workflow_results = [(w, w.result()) for w in workflows]

for workflow, result in workflow_results:
    results[workflow.name] = result.affinity_score.probability_binary

print(results)
