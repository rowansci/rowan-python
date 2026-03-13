from stjames import Method, Molecule

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."

workflow = rowan.submit_basic_calculation_workflow(
    initial_molecule=Molecule.from_smiles("CC(=C)C=C"),
    method=Method.OMOL25_CONSERVING_S,
    tasks=["energy"],
    mode="auto",
    engine="omol25",
    name="Isoprene Energy",
)

print(f"View workflow privately at: https://labs.rowansci.com/calculation/{workflow.uuid}")
result = workflow.result()
print(result)
