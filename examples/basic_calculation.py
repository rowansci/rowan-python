import time

from stjames import Method, Molecule

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

workflow = rowan.submit_basic_calculation_workflow(
    initial_molecule=Molecule.from_smiles("CC(=C)C=C"),
    method=Method.OMOL25_CONSERVING_S,
    tasks=["optimize"],
    mode="auto",
    engine="omol25",
    name="Isoprene Optimization",
    folder=folder,
)

print(f"View workflow privately at: https://labs.rowansci.com/calculation/{workflow.uuid}")

# Poll the calculation endpoint to watch optimization steps as they complete.
calc_uuid = None
while not workflow.done():
    partial = workflow.result(wait=False)
    if calc_uuid is None:
        calc_uuid = partial.calculation_uuid
    if calc_uuid:
        mols = rowan.retrieve_calculation_molecules(calc_uuid)
        print(f"  {len(mols)} opt steps, energy={mols[-1].get('energy') if mols else None}")
    time.sleep(3)

result = workflow.result()
print(result)
