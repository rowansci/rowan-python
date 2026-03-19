"""
Calculate hydrogen bond acceptor/donor strength using the Rowan API.

This workflow predicts pKBHX values of hydrogen-bond acceptors and pKa values
of hydrogen-bond donors using neural network potentials and r2SCAN-3c.

See documentation at: https://docs.rowansci.com/science/workflows/hydrogen-bond-basicity
"""

from stjames import Molecule

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

workflow = rowan.submit_hydrogen_bond_donor_acceptor_strength_workflow(
    initial_molecule=Molecule.from_smiles("CC(=O)N(C)C"),  # dimethylacetamide
    do_csearch=True,  # run conformer search (default)
    do_optimization=True,  # optimize structures (default)
    name="DMA hydrogen bond basicity",
    folder=folder,
)

print(
    f"View workflow privately at: https://labs.rowansci.com/hydrogen-bond-basicity/{workflow.uuid}"
)
result = workflow.result()
print(result)
# e.g. <HydrogenBondDonorAcceptorStrengthResult acceptors=2 donors=1>
