from pathlib import Path

import rowan

protein = rowan.upload_protein("1IEP receptor", Path("examples/data/1iep_receptorH.pdb"))

md_workflow = rowan.submit_protein_md_workflow(
    protein=protein,
    num_trajectories=1,
    simulation_time_ns=1,
    name="MD on 1IEP (apo)",
)

print(f"View MD workflow privately at: https://labs.rowansci.com/protein-md/{md_workflow.uuid}")
