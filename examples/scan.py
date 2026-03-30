"""
Run an scan calculation on a molecule using Rowan.

See documentation at: https://docs.rowansci.com/science/workflows/scan
"""

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

workflow = rowan.submit_scan_workflow(
    initial_molecule=rowan.Molecule.from_smiles("O"),
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
    folder=folder,
)

print(f"View workflow privately at: https://labs.rowansci.com/scan/{workflow.uuid}")
result = workflow.result()
print(result)
# e.g. <ScanResult points=5>
