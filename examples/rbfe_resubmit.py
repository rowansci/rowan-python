"""Extend a converged RBFE graph with an additional ligand via reuse.

A completed RBFE graph holds converged per-edge relative free energies (ddG) on
each ligand-pair edge. Passing it back as `seed_graph` freezes those edges and
their results, constructs only the edges incident to the new ligand (atom-mapping
plus alchemical-path setup), and restricts the subsequent FEP to those edges -- the
existing complex/solvent legs are not resimulated. A per-ligand relative binding
free energy for the new ligand (relative to the reference ligand the graph is
anchored on) is then recovered by re-solving the maximum-likelihood estimate over
the augmented edge set.
"""

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."

# UUID of a COMPLETED relative_binding_free_energy_perturbation workflow to extend.
FINISHED_RBFE_UUID = "your-completed-rbfe-uuid"

folder = rowan.get_folder("examples")

finished = rowan.retrieve_workflow(FINISHED_RBFE_UUID).result()
print(f"Finished study: {len(finished.ligands)} ligands, {len(finished.edges)} edges with ddG")

# Combine the finished study's ligands with the newly designed one(s) to score.
new_ligands = rowan.load_named_ligands("path/to/new_ligands.sdf")
ligands = {**finished.ligands, **new_ligands}

extended_graph = rowan.submit_relative_binding_free_energy_graph_workflow(
    ligands=ligands,
    seed_graph=finished.graph,
    folder=folder,
    name="RBFE Graph (extended)",
).result()

resubmit = rowan.submit_relative_binding_free_energy_perturbation_workflow(
    graph_result=extended_graph,
    protein=finished.protein,
    folder=folder,
    name="RBFE Perturbation (resubmit)",
)
print(
    "Resubmitted! View at: "
    f"https://labs.rowansci.com/relative-binding-free-energy-perturbation/{resubmit.uuid}"
)

result = resubmit.result()
for name in new_ligands:
    dg = result.ligand_dg_results[name]
    print(f"{name}: dG = {dg.dg:.2f} +/- {dg.dg_err:.2f} kcal/mol")
