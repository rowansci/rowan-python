from pathlib import Path

import stjames

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."

folder = rowan.get_folder("examples")

citalopram_analogues = [
    "CN(C)CCC[C@@]1(c2ccccc2)OCc2cc(C#N)ccc21",
    "CN(C)CCC[C@@]1(c2ccc(F)cc2)OCc2c(CC)c(C#N)ccc21",
    "CN(C)CCC[C@@]1(c2ccc(CCC)cc2)OCc2cc(C#N)ccc21",
]

data_dir = Path(__file__).parent / "data"
bound_pose = stjames.Molecule.from_file(str(data_dir / "citalopram_1iep.xyz"))

protein = rowan.upload_protein("1IEP receptor", data_dir / "1iep_receptorH.pdb")

workflow = rowan.submit_analogue_docking_workflow(
    analogues=citalopram_analogues,
    protein=protein,
    initial_molecule=bound_pose,
    folder=folder,
)

print(f"View MD workflow privately at: https://labs.rowansci.com/analogue-docking/{workflow.uuid}")
result = workflow.result()

print(result)
# e.g. <AnalogueDockingResult analogues=3 best=(-8.30, 'CN(C)CCC...')>
