import cctk
import rowan

rowan.api_key = "rowan-sk..."
client = rowan.Client()

# load molecule by name (cctk can also load in a variety of file formats)
molecule = cctk.Molecule.new_from_name("cyclobutane")

# run calculation remotely and return result
result = client.compute(
    molecule,
    name="opt cyclobutane",
    method="b97-3c",
    tasks=["optimize", "charge", "dipole"]
)

print(result)
