from stjames import Molecule

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."

workflow = rowan.submit_ion_mobility_workflow(
    Molecule.from_smiles("c1ccccn1"),
    name="pyridinium CCS",
    protonate=True,
)

print(f"View workflow privately at: https://labs.rowansci.com/ion_mobility/{workflow.uuid}")
result = workflow.result()

print(result.average_ccs)  # Å**2
