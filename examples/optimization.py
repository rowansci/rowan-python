"""
Run an optimization calculation on a molecule using Rowan.

See documentation at: https://docs.rowansci.com/science/quantum-chemistry/geometry-optimization
"""

import json

from stjames import Molecule

import rowan

# rowan.api_key = ""

# Run calculation remotely
result = rowan.compute(
    Molecule.from_smiles("O"),
    workflow_type="basic_calculation",
    name="Water Optimization",
    settings={
        "method": "GFN2-xTB",
        "tasks": ["optimize"],
    },
    engine="xtb",
)

print(json.dumps(result, indent=4))
