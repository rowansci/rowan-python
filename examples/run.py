import rowan
import stjames 

rowan.api_key = "rowan-sk..."

# load molecule by smiles (stjames can also load in a variety of file formats)
molecule = stjames.Molecule.from_smiles("C1CCC1") # cyclobutane

# run calculation remotely and return result
result = rowan.compute(
    molecule,
    name="opt cyclobutane",
    method="b97-3c",
    tasks=["optimize", "charge", "dipole"]
)

print(result)
