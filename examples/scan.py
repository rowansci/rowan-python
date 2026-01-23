"""
Run an scan calculation on a molecule using Rowan.

See documentation at: https://docs.rowansci.com/science/workflows/scan
"""

from stjames import Molecule

import rowan

# Set ROWAN_API_KEY environment variable to your API key or set rowan.api_key directly
# rowan.api_key = "rowan-sk..."

workflow = rowan.submit_scan_workflow(
    initial_molecule=Molecule.from_smiles("O"),
    name="Water Angle scan",
    scan_settings={
        "type": "angle",
        "atoms": [2, 1, 3],  # 1-indexed
        "start": 100,
        "stop": 110,
        "num": 5,
    },
    calculation_method="GFN2-xTB",
    calculation_engine="xtb",
)

print(f"View workflow privately at: https://labs.rowansci.com/workflow/{workflow.uuid}")
workflow.wait_for_result().fetch_latest(in_place=True)
print(workflow)
