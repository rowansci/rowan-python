"""
Calculate Bond-Dissociation Energies (BDE) with the Rowan API.

`mode` is a BDE method string (the level of theory):
- "omol25_conserving_s": neural network potential (default)
- "g_xtb//gfn2_xtb": semiempirical
- "r2scan3c//gfn2_xtb": DFT single point

Choose which bonds to break:
- all_CH=True / all_CX=True: every C-H / every C-X (carbon-halogen) bond
- fragment_indices=[[...]]: 1-indexed atoms of a fragment to split off (each fragment must
  connect to the rest by a single bond) -- e.g. [[9]] splits off one hydrogen,
  [[3, 9]] splits off a whole group like -OH

Use the bond finders to get indices: rowan.find_ch_bonds(mol), rowan.find_cx_bonds(mol),
or rowan.find_bonds(mol, element_a, element_b, distance_max).

See documentation at: https://docs.rowansci.com/science/workflows/bond-dissociation-energy
"""

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

mol = rowan.Molecule.from_smiles("CCCC")  # butane

# Break every C-H bond and report each one's BDE.
workflow = rowan.submit_bde_workflow(
    initial_molecule=mol,
    mode="omol25_conserving_s",
    all_CH=True,
    name="Butane C-H BDEs",
    folder=folder,
)

print(f"View workflow privately at: https://labs.rowansci.com/bde/{workflow.uuid}")
result = workflow.result()
print(result)
print(result.bdes)
# e.g. <BDEResult energy=-... Ha bdes=...>


# To target a specific bond instead of a whole class, pass its fragment atoms (1-indexed).
# This reaches bonds the all_CH / all_CX flags don't cover -- here, the ethanol O-H bond.
ethanol = rowan.Molecule.from_smiles("CCO")
oh_bonds = rowan.find_bonds(ethanol, element_a=8, element_b=1, distance_max=1.1)
oh_hydrogen = oh_bonds[0][1]  # the hydroxyl H

oh_workflow = rowan.submit_bde_workflow(
    initial_molecule=ethanol,
    mode="omol25_conserving_s",
    fragment_indices=[[oh_hydrogen]],  # split off that H -> the ethanol O-H BDE
    name="Ethanol O-H BDE",
    folder=folder,
)

print(f"View workflow privately at: https://labs.rowansci.com/bde/{oh_workflow.uuid}")
oh_result = oh_workflow.result()
print(oh_result)
print(oh_result.bdes)
