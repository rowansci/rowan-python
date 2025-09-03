from stjames import Molecule

import rowan

# Set ROWAN_API_KEY environment variable to your API key or set rowan.api_key directly
# rowan.api_key = "rowan-sk..."

result = rowan.submit_tautomer_search_workflow(
    initial_molecule=Molecule.from_smiles("O=C1C=CC=CN1"),
    mode="reckless",
    name="2-Pyridone Tautomers",
)


print(result.wait_for_result().fetch_latest(in_place=True))
