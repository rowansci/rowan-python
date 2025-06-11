"""
Run an scan calculation on a molecule using Rowan.

See documentation at: https://docs.rowansci.com/science/workflows/scan
"""

import json

from stjames import Molecule

import rowan

# rowan.api_key = ""

# Run calculation remotely
result = rowan.compute(
    Molecule.from_smiles("O"),
    workflow_type="scan",
    name="Water Angle scan",
    scan_settings={
        "type": "angle",
        "atoms": [2, 1, 3],  # 1-indexed
        "start": 100,
        "stop": 110,
        "num": 5,
    },
    calc_settings={
        "method": "GFN2-xTB",
        "tasks": ["optimize"],
    },
    calc_engine="xtb",
)

print(json.dumps(result, indent=4))
