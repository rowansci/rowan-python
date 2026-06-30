# Binding Affinity

## Input

Two submission modes — choose based on whether the ligand is already in the protein PDB.

**Mode 1 — holo protein:** the protein PDB already contains the bound ligand as a residue. Pass `ligand_residue_name` to identify which residue is the ligand; everything else is treated as receptor.

**Mode 2 — apo protein + external poses:** the protein has no bound ligand. Pass `ligand_structures` with one or more molecules that are already in the protein's coordinate frame (e.g. from a prior docking run). When scoring multiple ligands, this mode is preferred over submitting separate workflows — all poses share the same pocket geometry.

- Protein: a `rowan.Protein` or its UUID. Get one from the PDB with `rowan.create_protein_from_pdb_id(pdb_code)`, or upload your own with `rowan.upload_protein(name, path)`.
- `ligand_residue_name`: residue name of the bound ligand in a holo protein PDB (mode 1).
- `ligand_structures`: list of `StructureInput` molecules with 3D coordinates, already aligned to the protein (mode 2). Load from an SDF with `rowan.load_named_ligands(path)`.

Binding affinity is computed using SQM energies (PM6-D3H4X: COSMO geometry optimization followed by COSMO2 single-point, in water by default) on a truncated pocket around the ligand. When scoring multiple ligands, submitting them together in one workflow is preferred over separate workflows — all poses share the same pocket environment, making scores directly comparable.

## Example

```python
from pathlib import Path
import rowan

folder = rowan.get_folder("examples")
data_dir = Path("examples/data")

protein = rowan.upload_protein("TYK2", data_dir / "tyk2_structure.pdb")
ligands = rowan.load_named_ligands(data_dir / "tyk2_ligands.sdf")

workflow = rowan.submit_binding_affinity_workflow(
    protein=protein,
    ligand_structures=list(ligands.values()),
    name="Binding Affinity — TYK2 ligands",
    folder=folder,
)

result = workflow.result()
for name, score in zip(ligands.keys(), result.scores):
    print(f"{name}: {score.binding_affinity:.2f} kcal/mol (strain: {score.strain})")
```

## Settings

`binding_affinity_settings` accepts a `rowan.SinglePointEnergySettings` object:

- `multistage_opt_settings` (default PM6-D3H4X/COSMO optimization + PM6-D3H4X/COSMO2 single-point in water): a `rowan.MultiStageOptSettings` controlling ligand geometry optimization and energy evaluation.
- `truncation_radius` (default `6.0` Å): protein residues beyond this distance from the ligand are excluded from the calculation.

```python
rowan.submit_binding_affinity_workflow(
    protein=protein,
    ligand_residue_name="LIG",
    binding_affinity_settings=rowan.SinglePointEnergySettings(truncation_radius=8.0),
)
```

## Results

`result.scores` is a list of `BindingAffinityScore` objects in the same order as the input ligands:

- `binding_affinity`: binding affinity in kcal/mol (ΔE = E(complex) − E(protein_region) − E(ligand)).
- `strain`: energy difference between the input pose and the SQM-optimized pose, in kcal/mol.
