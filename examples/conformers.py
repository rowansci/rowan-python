"""
Calculate the conformers of a molecule using the Rowan API.

| Mode                                          | Reckless | Rapid | Careful | Meticulous |
|-----------------------------------------------|----------|-------|---------|------------|
| Conformer generation program                  | RDKit    | RDKit | CREST   | CREST      |
| Level of theory used for conformer generation | MMFF     | MMFF  | GFN-FF  | GFN2-xTB   |
| CREST Mode                                    |          |       | quick   | normal     |
| ETKDG number of initial conformers            | 200      | 300   |         |            |
| Initial energy cutoff (kcal/mol)              | 10       | 15    | 10      | 15         |
| RMSD similarity cutoff (Å)                    | 0.25     | 0.10  | 0.00    | 0.00       |
| xTB screening max number of conformers        | 50       | 100   | 150     | 500        |
| Final energy cutoff (kcal/mol)                | 5        | 5     | 5       | 10         |

Rapid is recommended for most work.

See documentation at: https://docs.rowansci.com/science/workflows/conformers
"""

import json

from stjames import Molecule

import rowan

# rowan.api_key = ""

# Run calculation remotely
result = rowan.compute(
    Molecule.from_smiles("CCOCC"),
    workflow_type="conformers",
    name="Diethyl ether conformers",
    mode="reckless",
)

print(json.dumps(result, indent=4))
