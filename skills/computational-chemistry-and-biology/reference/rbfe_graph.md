# RBFE graph

## Input

A named set of ligands with 3D coordinates, each a `rowan.Molecule`, `stjames.Molecule`, or RDKit `Mol` (`StructureInput`).

- `ligands`: a dict mapping ligand name to its structure. Usually load a congeneric series from an SDF, MOL, or MOL2 file with `rowan.load_named_ligands(path)`, which reads the name from each record's title field and preserves the 3D coordinates the graph builder uses. Giving each ligand an explicit SMILES is recommended (not required): the workflow maps atoms between ligands, and a SMILES improves the stability of that mapping rather than relying on connectivity inferred from the 3D coordinates. The mapping defines how each ligand pair is morphed during the perturbation, so its quality is crucial to the reliability of the downstream FEP free-energy results. Loading from an SDF that stores a SMILES per record is one easy way to get this.

This workflow turns the ligand list into an RBFE graph whose edges are the ligand pairs the relative binding free energy perturbation (FEP) workflow transforms between.

All ligands in a graph must share the same formal charge — the workflow raises a `ValueError` otherwise. RBFE across a charge boundary samples poorly, so run neutral and charged ligands as separate graphs (see the relative binding free energy perturbation reference for the rationale).

## Example

```python
from pathlib import Path
import rowan

folder = rowan.get_folder("examples")
data_dir = Path("examples/data")

# 3D coordinates from the SDF are used by the graph builder.
ligands = rowan.load_named_ligands(data_dir / "tyk2_ligands.sdf")

wf = rowan.submit_relative_binding_free_energy_graph_workflow(
    ligands=ligands,
    mode="greedy",
    folder=folder,
)

result = wf.result()
print(result)   # e.g. <RelativeBindingFreeEnergyGraphResult edges=18 ligands=16>

for edge in result.edges:
    print(f"  {edge.ligand_a} -> {edge.ligand_b}  score={edge.score:.3f}")
```

The resulting `RelativeBindingFreeEnergyGraphResult` is passed as `graph_result` to the relative binding free energy perturbation workflow.

## Extending a finished study with new ligands (resubmit)

The main use of `seed_graph` is to add ligands to a study that has **already run the full FEP**, without recomputing the existing edges. Pass the finished graph as `seed_graph` along with the full ligand set (old plus new); its edges are preserved verbatim — including the ddG already computed for them — and only edges connecting the new ligands are built. When you then resubmit the perturbation, only the new edges run.

```python
# Retrieve a completed perturbation; its graph already carries ddG on every edge.
finished = rowan.retrieve_workflow(FINISHED_RBFE_UUID).result()
ligands = finished.ligands
ligands["new-analog"] = new_ligand   # a rowan.Molecule / stjames.Molecule / RDKit Mol

# Rebuild the graph, seeding from the finished graph: existing edges (and ddG) kept,
# only edges for "new-analog" are added.
graph_result = rowan.submit_relative_binding_free_energy_graph_workflow(
    ligands=ligands,
    seed_graph=finished.graph,
    folder=folder,
).result()

# Resubmit the FEP: preserved edges are reused, only the new edges run.
rowan.submit_relative_binding_free_energy_perturbation_workflow(
    graph_result=graph_result,
    protein=finished.protein,  # reuse the same target the study was run against
    folder=folder,
)
```

`seed_graph` takes the RBFE graph from a completed result (`result.graph`). Every ligand referenced by a seed-graph edge must be present in `ligands`. See `examples/rbfe_resubmit.py` for the full flow.

## Settings

- `mode` (default `greedy`): graph construction strategy, `greedy` or `star_map`.
- `hub_compound_id` (default none): the ligand name used as the hub when `mode="star_map"`.
- `greedy_scoring` (default `best`): how greedy mode scores candidate edges. `best` tries both other methods and keeps the graph with the fewest dummy atoms; `jaccard` scores by Jaccard distance of the atom-mapping overlap; `dummy_atoms` scores by the number of dummy atoms.
- `greedy_k_min_cut` (default `3`): minimum number of edges per ligand node, i.e. how many edges can be cut before the graph disconnects. Higher gives a more redundant, robust graph at the cost of more FEP edges to run. Must be greater than 0.
- `refine_cutoff` (default none): optional maximum-common-substructure similarity cutoff for graph refinement.
- `seed_graph` (default none): the RBFE graph from a prior run (`result.graph`) to extend. Its edges (and their computed results) are preserved; only edges for newly added ligands are built. See "Extending a graph with new ligands" above.

## Result fields

- `edges`: the graph edges (the ligand pairs FEP will transform between); each has `ligand_a`, `ligand_b`, and `score`.
- `graph`: the constructed `RBFEGraph`; pass back as `seed_graph` to extend the study.
- `ligands`: the input ligand molecules, keyed by identifier.
