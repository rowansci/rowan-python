from stjames import Molecule

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
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

print(f"View workflow privately at: https://labs.rowansci.com/calculation/{workflow.uuid}")
result = workflow.result()
print(result)
