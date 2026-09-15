# Binding Affinity

## Input

Choose one input mode:

**Holo protein:** pass a `protein` containing the bound ligand and use `ligand_residue_name` to identify it.

**Protein + external poses:** pass a `protein` and one or more aligned `ligand_structures`. For ML scoring, include each molecule's `smiles` when known so bond orders can be assigned reliably.

**NESSO sequence/SMILES:** pass `protein_sequences` and `ligand_smiles` with `rowan.NessoAffinitySettings()`; no PDB is required.

NESSO also accepts the first two modes, but uses only the protein sequence from a PDB rather than its 3D coordinates.

- Protein: a `rowan.Protein` or its UUID. Create one with `rowan.create_protein_from_pdb_id(...)` or `rowan.upload_protein(...)`.
- `ligand_structures`: 3D molecules already aligned to the protein, commonly loaded with `rowan.load_named_ligands(...)`.
- `protein_sequences` and `ligand_smiles`: PDB-free inputs supported only by NESSO.

Four scoring methods, selected via `binding_affinity_settings`:

- `rowan.SinglePointEnergySettings()` (default): SQM energies (PM6-D3H4X: COSMO geometry optimization followed by COSMO2 single-point, in water by default) on a truncated pocket around the ligand. Result in kcal/mol. Works in modes 1/2 only.
- `rowan.GninaAffinitySettings()`: GNINA CNN affinity scoring. Result in log10(M). Works in modes 1/2.
- `rowan.AEVPLIGAffinitySettings()`: AEV-PLIG affinity scoring. Result in log10(M). Works in modes 1/2.
- `rowan.NessoAffinitySettings()`: NESSO affinity prediction. Result in log10(M). Works in modes 1/2/3 — the only method that supports mode 3.

## Example

```python
from pathlib import Path
import rowan

folder = rowan.get_folder("examples")
data_dir = Path("examples/data")

protein = rowan.upload_protein("TYK2", data_dir / "tyk2_structure.pdb")
ligands = rowan.load_named_ligands(data_dir / "tyk2_ligands.sdf")

workflow = rowan.submit_binding_affinity_workflow(
    protein=protein.uuid,
    ligand_structures=list(ligands.values()),
    name="Binding Affinity — TYK2 ligands",
    folder=folder,
)

result = workflow.result()
for name, score in zip(ligands.keys(), result.scores):
    if score is None:
        print(f"{name}: scoring failed")
        continue
    print(f"{name}: {score.binding_affinity:.2f} kcal/mol")
```

## Settings

`binding_affinity_settings` accepts one of four settings objects:

- `rowan.SinglePointEnergySettings` (default):
  - `multistage_opt_settings` (default PM6-D3H4X/COSMO optimization + PM6-D3H4X/COSMO2 single-point in water): a `rowan.MultiStageOptSettings` controlling ligand geometry optimization and energy evaluation.
  - `truncation_radius` (default `6.0` Å): protein residues beyond this distance from the ligand are excluded from the calculation.
- `rowan.GninaAffinitySettings`, `rowan.AEVPLIGAffinitySettings`, `rowan.NessoAffinitySettings`: no tunable parameters.

```python
rowan.submit_binding_affinity_workflow(
    protein=protein.uuid,
    ligand_residue_name="LIG",
    binding_affinity_settings=rowan.SinglePointEnergySettings(truncation_radius=8.0),
)

# GNINA / AEV-PLIG: same submission shape, no PDB-free mode.
rowan.submit_binding_affinity_workflow(
    protein=protein.uuid,
    ligand_structures=list(ligands.values()),
    binding_affinity_settings=rowan.GninaAffinitySettings(),
)

# NESSO: sequence/SMILES input, no PDB required.
rowan.submit_binding_affinity_workflow(
    protein_sequences=["ACDEFGHIK"],
    ligand_smiles=["CCO"],
    binding_affinity_settings=rowan.NessoAffinitySettings(),
)
```

## Results

`result.scores` is a list in the same order as the input ligands. Each entry is a
`BindingAffinityScore`, or `None` when that input failed:

- `binding_affinity`: binding affinity in kcal/mol (ΔE = E(complex) − E(protein_region) − E(ligand)) for `SinglePointEnergySettings`, or in log10(M) for `GninaAffinitySettings`, `AEVPLIGAffinitySettings`, and `NessoAffinitySettings`.
