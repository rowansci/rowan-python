"""Run crambin MD with GBn2 implicit water."""

import rowan

protein = rowan.create_protein_from_pdb_id("1CRN", name="crambin")
preparation = rowan.submit_protein_preparation_workflow(
    protein=protein.uuid, add_missing_method="pdbfixer"
)
prepared_protein_uuid = preparation.result().prepared_protein_uuid
if prepared_protein_uuid is None:
    raise ValueError("Protein preparation returned no prepared protein")

workflow = rowan.submit_protein_md_workflow(
    protein=prepared_protein_uuid,
    implicit_solvent="gbn2",
    pressure_atm=None,
    nonbonded_cutoff=None,
    water_buffer=None,
    ionic_strength_M=0.15,
    num_trajectories=1,
    equilibration_time_ns=0.1,
    simulation_time_ns=1,
    name="Crambin in implicit water",
)
print(f"View MD workflow: https://labs.rowansci.com/protein-md/{workflow.uuid}")
result = workflow.result()
print(f"Protein RMSD: {result.trajectories[0].protein_rmsd}")
