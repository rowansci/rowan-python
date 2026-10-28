# Pocket detection

## Input

A `rowan.Protein` (or its UUID string). Proteins are not SMILES. Load one and prepare it before submitting:

- From a PDB ID: `rowan.create_protein_from_pdb_id("1OTP", name="thymidine phosphorylase")`.
- From a local PDB file: `rowan.upload_protein("my protein", "path/to/file.pdb")`.

Pocket detection accepts any stored `rowan.Protein` or protein UUID. Protein preparation is recommended before submission; when chaining from it, passing `prepared_protein_uuid` avoids an unnecessary structure fetch.

Detects potential binding pockets on the protein structure using Pocketeer (default), or with mixed-solvent molecular dynamics (MixMD), which finds regions where small probe molecules accumulate across independent replicas.

## Example

```python
import rowan

folder = rowan.get_folder("examples")

protein = rowan.create_protein_from_pdb_id("1OTP", name="thymidine phosphorylase")
preparation_workflow = rowan.submit_protein_preparation_workflow(
    protein=protein.uuid,
    folder=folder,
)
prepared_protein_uuid = preparation_workflow.result().prepared_protein_uuid

wf = rowan.submit_pocket_detection_workflow(
    protein=prepared_protein_uuid,
    folder=folder,
)

result = wf.result()
print(f"Detected {len(result.pockets)} pocket(s)")
for i, pocket in enumerate(result.pockets):
    print(i, "score=", pocket.score, "volume=", pocket.volume)
    print("  center=", pocket.pocket_center, "residues=", pocket.residue_numbers)
```

## Settings

- `merge_distance` (default `1.75`): Pocketeer distance in angstrom for merging detected pocket spheres into a single pocket.
- `settings`: `rowan.PocketeerSettings(merge_distance=...)` (default) or `rowan.MixMDSettings(...)`. MixMD settings include `probe` (`rowan.Solvent.ACETONITRILE` or `rowan.Solvent.ISOPROPANOL`), `num_trajectories` (default `3`), `equilibration_time_ns` (default `5`), `simulation_time_ns` (default `15`), `enrichment_threshold` (default `5`, probe density relative to bulk), and the same force-field and simulation settings as protein MD. Set `merge_distance` on `PocketeerSettings` when passing `settings`.
- `validate_forcefield` (default `True`): for MixMD, check that the protein can be parameterized before submitting.

## Result fields

`result.pockets` is a list of detected pockets, each with:

- `score`: method-specific ranking score (larger is better); scores are not comparable between Pocketeer and MixMD.
- `volume`: pocket volume in cubic angstrom. For MixMD this is the volume of the enriched probe-density region, not the cavity.
- `pocket_center`: center of the axis-aligned bounding box, in angstrom.
- `pocket_sides`: side lengths of the bounding box, in angstrom.
- `residue_numbers`: residue numbers lining the pocket.
- `sphere_centers`, `sphere_radii`: centers and radii of the detected spheres, in angstrom.
- `mixmd`: for MixMD pockets, `density_center`, `peak_enrichment`, and `integrated_excess_density`; otherwise `None`.

For MixMD runs, `result.mixmd` lists each replica (`seed`, `num_frames`, `num_probes`, `largest_probe_cluster_fraction`, and artifact UUIDs), and `result.get_reference_protein()` returns the protein that defines pocket coordinates. Reliability warnings, such as fewer than three replicas or disagreement between replicas, appear in `result.messages`.

To dock into a detected pocket, pass its `pocket_center` and `pocket_sides` as the docking workflow's `pocket`: `pocket=[list(p.pocket_center), list(p.pocket_sides)]`.
