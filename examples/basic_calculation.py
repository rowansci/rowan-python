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

# Stream optimization steps as they complete; final iteration is the complete result.
for result in workflow.stream_result(poll_interval=3):
    if result.calculation_uuid:
        mols = rowan.retrieve_calculation_molecules(result.calculation_uuid)
        print(f"  {len(mols)} opt steps, energy={mols[-1].get('energy') if mols else None}")

print(result)
# e.g. <BasicCalculationResult energy=-76.234567 H>
