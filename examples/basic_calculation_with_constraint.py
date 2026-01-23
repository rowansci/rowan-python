from stjames import Molecule

import rowan

# Set ROWAN_API_KEY environment variable to your API key or set rowan.api_key directly
# rowan.api_key = "rowan-sk..."

workflow = rowan.submit_workflow(
    initial_molecule=Molecule.from_smiles("CCCC"),
    workflow_type="basic_calculation",
    name="Constrained Butane",
    workflow_data={
        "settings": {
            "method": "gfn2_xtb",
            "tasks": ["optimize"],
            "mode": "auto",
            "opt_settings": {
                "constraints": [
                    {
                        "atoms": [4, 3, 2, 1],
                        "constraint_type": "dihedral",
                        "value": 0,
                    },
                ]
            },
        },
    },
)

print(f"View workflow privately at: https://labs.rowansci.com/workflow/{workflow.uuid}")
workflow.wait_for_result().fetch_latest(in_place=True)
print(workflow)
