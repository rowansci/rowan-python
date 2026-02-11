"""
Calculate spin state energies using the Rowan API.

This workflow predicts the lowest energy spin state by running multistage
optimizations at different spin multiplicities.

Available modes:
- "reckless": Fastest, least accurate
- "rapid": Good balance of speed and accuracy (default)
- "careful": More accurate, slower
- "meticulous": Most accurate, slowest

See documentation at: https://docs.rowansci.com/science/workflows/spin-states
"""

from stjames import Molecule

import rowan

# Set ROWAN_API_KEY environment variable to your API key or set rowan.api_key directly
# rowan.api_key = "rowan-sk..."

# Methylene (CH2) - classic spin states example, triplet is ground state
workflow = rowan.submit_spin_states_workflow(
    initial_molecule=Molecule.from_smiles("[CH2]"),
    states=[1, 3],  # singlet vs triplet
    mode="rapid",
    name="Methylene spin states",
)

print(f"View workflow at: https://labs.rowansci.com/spin-states/{workflow.uuid}")
workflow.wait_for_result().fetch_latest(in_place=True)
print(workflow.data)
