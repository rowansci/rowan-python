from stjames import Molecule

import rowan

# Set ROWAN_API_KEY environment variable to your API key or set rowan.api_key directly
# rowan.api_key = "rowan-sk..."

workflow = rowan.submit_ion_mobility_workflow(
    Molecule.from_smiles("c1ccccn1"),
    name="pyridinium CCS",
    protonate=True,
)

print(f"View workflow privately at: https://labs.rowansci.com/workflow/{workflow.uuid}")
workflow.wait_for_result().fetch_latest(in_place=True)

print(workflow.data["average_ccs"])  # Å**2
