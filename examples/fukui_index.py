import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

workflow = rowan.submit_fukui_workflow(
    initial_molecule=rowan.Molecule.from_smiles("C1=CC=C(C=C1)C(=O)O"),
    optimization_method="gfn2_xtb",
    fukui_method="gfn1_xtb",
    name="Benzoic Acid Fukui",
    folder=folder,
)


print(f"View workflow privately at: https://labs.rowansci.com/fukui/{workflow.uuid}")
result = workflow.result()
print(result)
# e.g. <FukuiResult global_electrophilicity_index=2.34 eV>
