import stjames

import rowan

# Set ROWAN_API_KEY environment variable to your API key or set rowan.api_key directly
# rowan.api_key = "rowan-sk..."

workflow = rowan.submit_workflow(
    initial_molecule=stjames.Molecule.from_smiles("C1CCC1"),  # cyclobutane
    workflow_type="multistage_opt",
    name="Multistage optimization cyclobutane",
    workflow_data={
        "mode": "reckless",
    },
)

print(f"View workflow privately at: https://labs.rowansci.com/multistage-opt/{workflow.uuid}")
workflow.wait_for_result().fetch_latest(in_place=True)
print(workflow)
