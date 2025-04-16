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

See documentation at: https://docs.rowansci.com/science/workflows/bond-dissociation-energy
"""

import json

from stjames import Molecule

import rowan

# rowan.api_key = ""

# Run calculation remotely
result = rowan.compute(
    Molecule.from_smiles("CCCC"),
    workflow_type="bde",
    name="Butane BDE",
    mode="reckless",
    all_CH=True,
)

print(json.dumps(result, indent=4))
