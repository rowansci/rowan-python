from stjames import Molecule

import rowan

# Set ROWAN_API_KEY environment variable to your API key or set rowan.api_key directly
# rowan.api_key = "rowan-sk..."

result = rowan.submit_fukui_workflow(
    initial_molecule=Molecule.from_smiles("C1=CC=C(C=C1)C(=O)O"),
    optimization_method="gfn2_xtb",
    fukui_method="gfn1_xtb",
    name="Benzoic Acid Fukui",
)


print(result.wait_for_result().fetch_latest(in_place=True))
