from pathlib import Path

import stjames

import rowan

citalopram_analogues = [
    "CN(C)CCC[C@@]1(c2ccccc2)OCc2cc(C#N)ccc21",
    "CN(C)CCC[C@@]1(c2ccc(F)cc2)OCc2c(CC)c(C#N)ccc21",
    "CN(C)CCC[C@@]1(c2ccc(CCC)cc2)OCc2cc(C#N)ccc21",
]

bound_pose = stjames.Molecule.from_file("examples/data/citalopram_1iep.xyz")

protein = rowan.upload_protein("1IEP receptor", Path("examples/data/1iep_receptorH.pdb"))

workflow = rowan.submit_analogue_docking_workflow(
    analogues=citalopram_analogues,
    protein=protein,
    initial_molecule=bound_pose,
)

print(f"View MD workflow privately at: https://labs.rowansci.com/analogue-docking/{workflow.uuid}")
workflow.wait_for_result().fetch_latest(in_place=True)

# print ligand RMSD by frame
print(workflow.data)
