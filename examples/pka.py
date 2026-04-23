"""
Calculate the pKa of phenol using the Rowan API.

Experimental value ≈ 9.99

Four methods are available:
- aimnet2_wagen2024: AIMNet2-based; requires 3D structure; water only.
- gxtb_wagen2026: g-xTB-based; requires 3D structure; water only; full periodic table.
- chemprop_nevolianis2025: Chemprop-based; requires SMILES; several solvents supported.
- starling: SMILES-based; water only.

See documentation at: https://docs.rowansci.com/science/workflows/pka
"""

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples/pka")

# 3D structure-based methods — pass a Molecule with coordinates
mol = rowan.Molecule.from_smiles("c1ccccc1O")

workflow = rowan.submit_pka_workflow(
    initial_molecule=mol,
    method="aimnet2_wagen2024",
    name="Phenol pKa (aimnet2_wagen2024)",
    folder=folder,
)
print(f"View folder privately at: https://labs.rowansci.com/folder/{folder.uuid}")
print(f"View at: https://labs.rowansci.com/pka/{workflow.uuid}")
result = workflow.result()
print(result)
# e.g. <pKaResult acid=9.87>

workflow2 = rowan.submit_pka_workflow(
    initial_molecule=mol,
    method="gxtb_wagen2026",
    name="Phenol pKa (gxtb_wagen2026)",
    folder=folder,
)
print(f"View at: https://labs.rowansci.com/pka/{workflow2.uuid}")
result2 = workflow2.result()
print(result2)
# e.g. <pKaResult acid=9.92>

# SMILES-based methods — pass a SMILES string directly
smiles = "c1ccccc1O"

workflow3 = rowan.submit_pka_workflow(
    initial_molecule=smiles,
    method="chemprop_nevolianis2025",
    name="Phenol pKa (chemprop_nevolianis2025)",
    folder=folder,
)
print(f"View at: https://labs.rowansci.com/pka/{workflow3.uuid}")
result3 = workflow3.result()
print(result3)
# e.g. <pKaResult acid=9.75>

workflow4 = rowan.submit_pka_workflow(
    initial_molecule=smiles,
    method="starling",
    name="Phenol pKa (starling)",
    folder=folder,
)
print(f"View at: https://labs.rowansci.com/pka/{workflow4.uuid}")
result4 = workflow4.result()
print(result4)
# e.g. <pKaResult acid=9.81>
