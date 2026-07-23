# Protein binder design

## Input

A target to bind and a specification of the binder to design, passed as a `binder_design_input` dict in the BoltzGen format. It holds up to three entity lists; give the target through any of them, plus a designable binder:

- `protein_entities`: protein chains, each an `id` and a `sequence`. A sequence is either a fixed amino-acid string (a target chain) or a designable-region spec (the binder being designed).
- `ligand_entities`: small-molecule targets, each an `id` and a `smiles`.
- `file_entities`: 3D protein-structure targets, each referencing an uploaded protein by `uuid` (from `rowan.upload_protein` or `rowan.create_protein_from_pdb_id`), with optional region selections marking which residues to keep, exclude, or design.

You always need at least one `protein_entities` sequence with a designable region (the binder being designed), plus a target for it to bind (a fixed sequence, a structure, or a ligand). A target on its own, such as just a ligand SMILES, is not enough: without a designable protein sequence there is nothing to design, and submitting raises a `ValueError`.

Designable-region syntax for a binder sequence:

- `"140..180"`: design a protein of 140 to 180 residues.
- `"9A9C"`: design 9 residues, fix an alanine, design 9 more, fix a cysteine.
- `"MKTAYIAKQ"`: a fixed amino-acid sequence, with no designable region.

This workflow generates, filters, and ranks binders for a protein or small-molecule target.

## Example

```python
import rowan

folder = rowan.get_folder("examples")

# Design a protein to bind brilacidin, an antimicrobial small molecule.
binder_design_input = {
    "protein_entities": [
        {"id": "A", "sequence": "140..180"},  # design a 140-180 residue protein binder
    ],
    "ligand_entities": [
        {
            "id": "B",
            "smiles": "C1CNC[C@@H]1OC2=C(C=C(C=C2NC(=O)C3=CC(=NC=N3)C(=O)NC4=CC(=CC(=C4O[C@@H]5CCNC5)NC(=O)CCCCN=C(N)N)C(F)(F)F)C(F)(F)F)NC(=O)CCCCN=C(N)N",
        },
    ],
}

wf = rowan.submit_protein_binder_design_workflow(
    binder_design_input=binder_design_input,
    protocol="protein-small_molecule",
    folder=folder,
)

result = wf.result()
for binder in result.generated_binders:  # sorted by quality score
    print(binder.sequence, binder.iptm, binder.bound_structure_uuid)
```

## Settings

- `protocol` (default `protein-anything`): the design protocol, one of:
  - `protein-anything`: design a protein to bind a protein.
  - `peptide-anything`: design a peptide to bind a protein.
  - `protein-small_molecule`: design a protein to bind a small molecule.
  - `nanobody-anything`: design a single-domain antibody (nanobody) binder.
- `num_designs` (default `10`): number of designs to generate. The default is a small test; production runs typically use 10,000 to 60,000.
- `budget` (default `2`): how many designs are returned in the final diversity-optimized set, selected during BoltzGen's filtering step.

## Result fields

- `generated_binders`: the designs, sorted by quality score. Each has:
  - `sequence`: the designed binder's amino-acid sequence.
  - `bound_structure_uuid`: UUID of the predicted binder-target complex structure (feed it into a downstream workflow such as MD or docking).
  - `scores`: BoltzGen's metrics for the design (`iptm`, `design_ptm`, `quality_score`, `bb_rmsd`, and more). `iptm` and `quality_score` are also exposed directly on the binder for convenience.
