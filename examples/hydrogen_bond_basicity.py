"""
Calculate hydrogen bond acceptor/donor strength using the Rowan API.

This workflow predicts pKBHX values of hydrogen-bond acceptors and pKa values
of hydrogen-bond donors using neural network potentials and r2SCAN-3c.

See documentation at: https://docs.rowansci.com/science/workflows/hydrogen-bond-basicity
"""

from stjames import Molecule

import rowan

# Set ROWAN_API_KEY environment variable to your API key or set rowan.api_key directly
# rowan.api_key = "rowan-sk..."

workflow = rowan.submit_hydrogen_bond_basicity_workflow(
    initial_molecule=Molecule.from_smiles("CC(=O)N(C)C"),  # dimethylacetamide
    do_csearch=True,   # run conformer search (default)
    do_optimization=True,  # optimize structures (default)
    name="DMA hydrogen bond basicity",
)

print(f"View workflow at: https://labs.rowansci.com/hydrogen-bond-basicity/{workflow.uuid}")
workflow.wait_for_result().fetch_latest(in_place=True)
print(workflow.data)
