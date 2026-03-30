"""
Calculate the Bond-Dissociation Energy (BDE) of a molecule using the Rowan API.

| Calculation |   Initial    |    Final     |          Single        | Optimizes  |
| Mode        | Optimization | Optimization |          Point         | Fragments? |
|-------------|--------------|--------------|------------------------|------------|
| Reckless    |              |    GFN-FF    |        GFN2-xTB        |     No     |
| Rapid       |              |   GFN2-xTB   |        r²SCAN-3c       |     Yes    |
| Careful     |              |   r²SCAN-3c  |        ωB97X-3c        |     Yes    |
| Meticulous  |  r²SCAN-3c   |   ωB97X-3c   | ωB97M-D3BJ/def2-TZVPPD |     Yes    |

Rapid is recommended for most work.

Available modes:
- "reckless": Fastest, least accurate
- "rapid": Good balance of speed and accuracy (default)
- "careful": More accurate, slower
- "meticulous": Most accurate, slowest

See documentation at: https://docs.rowansci.com/science/workflows/bond-dissociation-energy
"""

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

workflow = rowan.submit_bde_workflow(
    initial_molecule=rowan.Molecule.from_smiles("CCCC"),
    mode="rapid",
    all_CH=True,
    name="Butane BDE",
    folder=folder,
)

print(f"View workflow privately at: https://labs.rowansci.com/bde/{workflow.uuid}")
result = workflow.result()
print(result)
# e.g. <BDEResult energy=-152.234567 Ha bdes=8>
