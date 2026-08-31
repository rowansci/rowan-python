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

# Cofolding predictions lack hydrogens — prepare and retain the ligand for MD
preparation_workflow = rowan.submit_protein_preparation_workflow(
    protein=cofolding_result.predicted_refined_structure_uuid,
    retain_non_polymer={"LIG": ligand},
    name="Prepare cofolded CDK2 complex",
    folder=folder,
)
prepared_protein_uuid = preparation_workflow.result().prepared_protein_uuid

md_workflow = rowan.submit_pose_analysis_md_workflow(
    protein=prepared_protein_uuid,
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

trajectory = md_result.trajectories[0]
print(f"Ligand RMSD: {trajectory.ligand_rmsd}")
print(f"Protein RMSD: {trajectory.protein_rmsd}")
print(f"Protein RMSF: {trajectory.rmsf}")
print(f"Potential energy: {trajectory.potential_energy} Hartree")
print(f"MM/GBSA scores: {trajectory.mmgbsa_scores} kcal/mol")
