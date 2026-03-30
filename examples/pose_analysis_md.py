import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

ligand = "CCC(C)(C)NC1=NCC2(CCC(=O)C2C)N1"

cofolding_workflow = rowan.submit_protein_cofolding_workflow(
    initial_protein_sequences=[
        "MENFQKVEKIGEGTYGVVYKARNKLTGEVVALKKIRLDTETEGVPSTAIREISLLKELNHPNIVKLLDVIHTENKLYLVFEFLHQDLKKFMDASALTGIPLPLIKSYLFQLLQGLAFCHSHRVLHRDLKPQNLLINTEGAIKLADFGLARAFGVPVRTYTHEVVTLWYRAPEILLGCKYYSTAVDIWSLGCIFAEMVTRRALFPGDSEIDQLFRIFRTLGTPDEVVWPGVTSMPDYKPSFPKWARQDFSKVVPPLDEDGRSLLSQMLHYDPNKRISAKAALAHPFFQDVTKPVPHLRL"
    ],
    initial_smiles_list=[ligand],
    ligand_binding_affinity_index=0,
    name=f"Cofolding {ligand}",
    do_pose_refinement=True,
    folder=folder,
)

print(
    f"View cofolding workflow privately at: https://labs.rowansci.com/protein-cofolding/{cofolding_workflow.uuid}"
)
cofolding_result = cofolding_workflow.result()

# Cofolding predictions lack hydrogens — prepare but keep the ligand for MD
protein = rowan.Protein(uuid=cofolding_result.predicted_refined_structure_uuid)
protein.prepare(remove_heterogens=False)

md_workflow = rowan.submit_pose_analysis_md_workflow(
    protein=protein,
    initial_smiles=ligand,
    num_trajectories=1,
    simulation_time_ns=1,
    name="Downstream molecular dynamics",
    folder=folder,
)

print(
    f"View MD workflow privately at: https://labs.rowansci.com/pose-analysis-md/{md_workflow.uuid}"
)
md_result = md_workflow.result()

# print ligand RMSD by frame
print(md_result.trajectories[0].ligand_rmsd)
