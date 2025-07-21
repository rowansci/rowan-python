import stjames

import rowan

# Set ROWAN_API_KEY environment variable to your API key or set rowan.api_key directly
# rowan.api_key = "rowan-sk..."

result = rowan.submit_workflow(
    initial_molecule=stjames.Molecule.from_smiles("C1CCC1"),  # cyclobutane
    workflow_type="multistage_opt",
    name="Multistage optimization cyclobutane",
    workflow_data={
        "mode": "reckless",
    },
)


print(result)
