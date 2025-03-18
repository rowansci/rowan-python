from pprint import pprint

import stjames

import rowan

# Set ROWAN_API_KEY environment variable to your API key or set rowan.api_key directly
# rowan.api_key = "rowan-sk..."

# load molecule by SMILES (stjames can load a variety of file formats)
molecule = stjames.Molecule.from_smiles("C1CCC1")  # cyclobutane

# run calculation remotely and return result
result = rowan.compute(
    molecule,
    workflow_type="multistage_opt",
    name="Multistage optimization cyclobutane",
    mode="reckless",
)

pprint(result)
