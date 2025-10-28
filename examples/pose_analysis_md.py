import rowan

ligand = "CCC(C)(C)NC1=NCC2(CCC(=O)C2C)N1"

cofolding_workflow = rowan.submit_protein_cofolding_workflow(
    initial_protein_sequences=[
        "MENFQKVEKIGEGTYGVVYKARNKLTGEVVALKKIRLDTETEGVPSTAIREISLLKELNHPNIVKLLDVIHTENKLYLVFEFLHQDLKKFMDASALTGIPLPLIKSYLFQLLQGLAFCHSHRVLHRDLKPQNLLINTEGAIKLADFGLARAFGVPVRTYTHEVVTLWYRAPEILLGCKYYSTAVDIWSLGCIFAEMVTRRALFPGDSEIDQLFRIFRTLGTPDEVVWPGVTSMPDYKPSFPKWARQDFSKVVPPLDEDGRSLLSQMLHYDPNKRISAKAALAHPFFQDVTKPVPHLRL"
    ],
    initial_smiles_list=[ligand],
    ligand_binding_affinity_index=0,
    name=f"Cofolding {ligand}",
    do_pose_refinement=True,
)

cofolding_workflow.wait_for_result().fetch_latest(in_place=True)

md_workflow = rowan.submit_pose_analysis_md_workflow(
    protein=cofolding_workflow.data["predicted_refined_structure_uuid"],
    initial_smiles=ligand,
    num_trajectories=1,
    simulation_time_ns=1,
    name="Downstream molecular dynamics",
)

md_workflow.wait_for_result().fetch_latest(in_place=True)

# print ligand RMSD by frame
print(md_workflow.data["trajectories"][0]["rmsd"])
