# Protein prep

Preparation is recommended before protein workflows, especially for experimental structures from the PDB. X-ray crystal structures commonly have missing residues and heavy atoms in unresolved regions, no hydrogens, and crystallization heterogens, all of which `prepare()` fixes. This file covers loading a protein, `prepare`, and `validate_protein_forcefield`.

## Loading a protein

Two ways to get a `rowan.Protein`:

- From a PDB code: `rowan.create_protein_from_pdb_id(code, name=None, project_uuid=None)` — name defaults to the PDB ID. Warns if the structure has multiple chains — crystal structures often have duplicate chains, so select the one you want with `protein.select_chains([...])` before running any workflow.
- From a local PDB file: `rowan.upload_protein(name, file_path, project_uuid=None)`

Both return a `rowan.Protein` you can then prepare.

## prepare()

The main entry point. Runs PDBFixer + OpenMM server-side and blocks until done.

```python
protein = rowan.create_protein_from_pdb_id("1IEP")

protein.prepare(
    find_missing_residues=True,  # identify and model missing residues
    add_missing_atoms=True,  # add missing heavy atoms
    remove_heterogens=True,  # remove ligands, salts, and other non-protein residues
    keep_waters=False,  # preserve waters when removing heterogens
    remove_hydrogens=False,  # remove existing hydrogens before re-adding
    remove_invalid_hydrogens=False,  # remove hydrogens not matching the forcefield template
    add_hydrogens=True,  # add missing hydrogens
    add_hydrogen_ph=7.0,  # pH used to set protonation states
    optimize_hydrogens=True,  # optimize hydrogen positions with OpenMM
    poll_interval=10.0,  # seconds between status checks
    timeout=300.0,  # seconds to wait before raising
)
```

Every argument above is shown at its default; `protein.prepare()` with no arguments does the same thing. The default prep fixes nonstandard residues, adds missing heavy atoms, removes heterogens, adds hydrogens at pH 7, and optimizes hydrogen positions, then blocks until the server reports done. (Nonstandard-residue replacement always runs; there is no flag for it.)

Some options only take effect alongside another:

- `keep_waters` only matters when `remove_heterogens=True`.
- `add_hydrogen_ph` and `optimize_hydrogens` only matter when `add_hydrogens=True`.
- `remove_hydrogens` removes all hydrogens, so also setting `remove_invalid_hydrogens` is redundant.

Common adjustments:

- **Keep a bound ligand** (for pose-analysis MD): `remove_heterogens=False`.
- **Keep structural waters**: `keep_waters=True`.
- **Different protonation pH**: `add_hydrogen_ph` (default `7.0`).
- **Input already has good hydrogens**: `add_hydrogens=False`.

Raises `RuntimeError` if preparation fails, is stopped, or times out.

## validate_protein_forcefield()

Server-side check that the protein can be parameterized by the MD forcefield. Call before any MD workflow (protein MD, pose-analysis MD, RBFE perturbation) to catch parameterization issues early.

```python
protein.validate_protein_forcefield(exclude_residue_names=None)
```

Ligand residues (`LIG`) are always excluded automatically. Pass `exclude_residue_names=["RES1", "RES2"]` for other residues to skip.

## Troubleshooting

If `validate_protein_forcefield()` fails, the usual fix is to re-prepare with `remove_invalid_hydrogens=True`, which strips hydrogens the validator did not like before re-adding them:

```python
protein.prepare(remove_invalid_hydrogens=True)
protein.validate_protein_forcefield()
```

## Pattern by use case

- **Docking, batch docking, analogue docking**: `prepare()` with defaults is fine.
- **Protein MD, pose-analysis MD, RBFE perturbation**: `prepare()` then `validate_protein_forcefield()` (for pose-analysis MD, use `prepare(remove_heterogens=False)` to keep the bound ligand for the MD run). On validation failure, re-`prepare(remove_invalid_hydrogens=True)` and validate again.
