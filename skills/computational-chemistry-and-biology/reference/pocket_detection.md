# Pocket detection

## Input

A `rowan.Protein` (or its UUID string). Proteins are not SMILES. Load one and prepare it before submitting:

- From a PDB ID: `rowan.create_protein_from_pdb_id("thymidine phosphorylase", "1OTP", project_uuid=rowan.default_project().uuid)`.
- From a local PDB file: `rowan.upload_protein("my protein", "path/to/file.pdb")`.

Then call `protein.prepare()`, which runs PDBFixer to fix nonstandard residues and add missing atoms and hydrogens, and waits for it to finish.

Detects potential binding pockets on the protein structure using Pocketeer.

## Example

```python
import rowan

folder = rowan.get_folder("examples")

protein = rowan.create_protein_from_pdb_id(
    "thymidine phosphorylase", "1OTP", project_uuid=rowan.default_project().uuid
)
protein.prepare()   # fix residues, add missing atoms/hydrogens; blocks until done

wf = rowan.submit_pocket_detection_workflow(
    protein=protein,
    folder=folder,
)

result = wf.result()
print(f"Detected {len(result.pockets)} pocket(s)")
for i, pocket in enumerate(result.pockets):
    print(i, "score=", pocket.score, "volume=", pocket.volume)
    print("  center=", pocket.pocket_center, "residues=", pocket.residue_numbers)
```

## Settings

- `merge_distance` (default `1.75`): distance in angstrom for merging detected pocket spheres into a single pocket.

## Result fields

`result.pockets` is a list of detected pockets, each with:

- `score`: druggability/quality score (larger is better).
- `volume`: pocket volume in cubic angstrom.
- `pocket_center`: center of the axis-aligned bounding box, in angstrom.
- `pocket_sides`: side lengths of the bounding box, in angstrom.
- `residue_numbers`: residue numbers lining the pocket.
- `sphere_centers`, `sphere_radii`: centers and radii of the detected spheres, in angstrom.

To dock into a detected pocket, pass its `pocket_center` and `pocket_sides` as the docking workflow's `pocket`: `pocket=[list(p.pocket_center), list(p.pocket_sides)]`.
