import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

workflow = rowan.submit_protein_cofolding_workflow(
    initial_protein_sequences=[
        "ASKGTSHEAGIVCRITKPALLVLNHETAKVIQTAFQRASYPDITGEKAMMLLGQVKYGLHNIQISHLSIASSQVELVEAKSIDVSIQDVSVVFKGTLKYGYTTAWWLGIDQSIDFEIDSAIDLQINTQLTADSGRVRTDAPDCYLSFHKLLLHLQGEREPGWIKQLFTNFISFTLKLVLKGQICKEINVISNIMADFVQTRAASILSDGDIGVDISLTGDPVITASYLESHHKGHFIYKDVSEDLPLPTFSPTLLGDSRMLYFWFSERVFHSLAKVAFQDGRLMLSLMGDEFKAVLETWGFNTNQEIFQEVVGGFPSQAQVTVHCLKMPKISCQNKGVVVDSSVMVKFLFPRPDQQHSVAYTFEEDIVTTVQASYSKKKLFLSLLDFQITPKTVSNLTESSSESIQSFLQSMITAVGIPEVMSRLEVVFTALMNSKGVSLFDIINPEIITRDGFLLLQMDFGFPEHLLVDFLQSLS"
    ],
    initial_smiles_list=[
        "CCOC(=O)N1c2ccc(C(F)(F)F)cc2[C@@H](N(Cc2cc(C(F)(F)F)cc(C(F)(F)F)c2)C(=O)OC)C[C@H]1CC"
    ],
    ligand_binding_affinity_index=0,
    name="Torcetrapib Cofolding",
    do_pose_refinement=True,
    compute_strain=False,
    folder=folder,
)

print(f"View workflow privately at: https://labs.rowansci.com/protein-cofolding/{workflow.uuid}")
result = workflow.result()
print(result)
