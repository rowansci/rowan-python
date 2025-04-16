"""
Calculate the pKa of the conjugate acid of pyridine using the Rowan API.

Experimental value ≈5.23

See documentiation at: https://docs.rowansci.com/science/workflows/pka
and preprint at: https://rowansci.com/publications/pka-prediction
"""

from stjames import Molecule

import rowan

# rowan.api_key = ""

# Run calculation remotely
result = rowan.compute(
    Molecule.from_smiles("n1ccccc1"),
    workflow_type="pka",
    name="Pyridine pKa",
    mode="reckless",
)

logfile = result["object_logfile"]
strongest_base = result["object_data"]["strongest_base"]

print(f"""\
{logfile}

pKa of conjugate acid: {strongest_base:.2f}""")
