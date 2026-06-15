# RBFE graph

## Input

A named set of ligands with 3D coordinates, each a `rowan.Molecule`, `stjames.Molecule`, or RDKit `Mol` (`StructureInput`).

- `ligands`: a dict mapping ligand name to its structure. Usually load a congeneric series from an SDF, MOL, or MOL2 file with `rowan.load_named_ligands(path)`, which reads the name from each record's title field and preserves the 3D coordinates the graph builder uses. Giving each ligand an explicit SMILES is recommended (not required): the workflow maps atoms between ligands, and a SMILES improves the stability of that mapping rather than relying on connectivity inferred from the 3D coordinates. The mapping defines how each ligand pair is morphed during the perturbation, so its quality is crucial to the reliability of the downstream FEP free-energy results. Loading from an SDF that stores a SMILES per record is one easy way to get this.

This workflow turns the ligand list into a perturbation graph whose edges are the ligand pairs the relative binding free energy perturbation (FEP) workflow transforms between.

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

## Settings

- `mode` (default `greedy`): graph construction strategy, `greedy` or `star_map`.
- `hub_compound_id` (default none): the ligand name used as the hub when `mode="star_map"`.
- `greedy_scoring` (default `best`): how greedy mode scores candidate edges. `best` tries both other methods and keeps the graph with the fewest dummy atoms; `jaccard` scores by Jaccard distance of the atom-mapping overlap; `dummy_atoms` scores by the number of dummy atoms.
- `greedy_k_min_cut` (default `3`): minimum number of edges per ligand node, i.e. how many edges can be cut before the graph disconnects. Higher gives a more redundant, robust graph at the cost of more FEP edges to run. Must be greater than 0.
- `refine_cutoff` (default none): optional maximum-common-substructure similarity cutoff for graph refinement.

## Result fields

- `edges`: the graph edges (the ligand pairs FEP will transform between); each has `ligand_a`, `ligand_b`, and `score`.
- `graph`: the constructed perturbation graph as a dict.
- `ligands`: the input ligand molecules, keyed by identifier.
