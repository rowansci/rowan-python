# ruff: noqa
from stjames import Molecule, Method

import rowan

# rowan.api_key = ""


result = rowan.submit_basic_calculation_workflow(
    initial_molecule=Molecule.from_smiles("CC(=C)C=C"),
    method=Method.OMOL25_CONSERVING_S,
    tasks=["energy"],
    mode="auto",
    engine="omol25",
    name="Isoprene Energy",
)

result.wait_for_result()
result.fetch_latest(in_place=True)

print(result)
