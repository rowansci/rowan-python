# Protein preparation

Protein preparation is recommended before protein workflows, especially for experimental structures from the PDB. X-ray crystal structures commonly have missing residues and heavy atoms in unresolved regions, no hydrogens, and crystallization heterogens. The protein-preparation workflow repairs the structure, assigns protonation states, and returns a new prepared protein without modifying the input.

## Loading a protein

Two ways to get a `rowan.Protein`:

- From a PDB code: `rowan.create_protein_from_pdb_id(code, name=None, project_uuid=None)` — name defaults to the PDB ID. Warns if the structure has multiple chains.
- From a local PDB file: `rowan.upload_protein(name, file_path, project_uuid=None)`

Both return a `rowan.Protein` you can then prepare.

## Selecting chains

Inspect a multi-chain structure before preparation. Crystal structures may contain duplicate copies of the target, unrelated peptides, antibodies, or other crystallization partners. Keep every chain needed for the intended binding site or biological assembly; do not reduce an interface or required oligomer to one chain. Chain order is not guaranteed, so select chains by their IDs rather than taking the first entry.

`select_chains()` returns a new protein record and leaves the original unchanged:

```python
protein = rowan.create_protein_from_pdb_id("4YHF")
print(protein.chains)

# Chain A contains BTK and its bound inhibitor in this structure.
protein = protein.select_chains(["A"])
```

Select chains before protein preparation. Non-polymer residues, waters, and branched entities assigned to discarded chains are removed with them. Check the PDB entry or inspect the structure when the correct chain is not obvious.

## Example

The defaults use Boltz-2 to add missing structure, cap termini with ACE/NME, protonate at pH 7.4 with OpenMM, retain existing protonation where possible, and keep sodium, chloride, and magnesium ions.

```python
import rowan

folder = rowan.get_folder("examples")
protein = rowan.create_protein_from_pdb_id("1CRN")

wf = rowan.submit_protein_preparation_workflow(
    protein=protein.uuid,
    folder=folder,
)
result = wf.result()
prepared_protein_uuid = result.prepared_protein_uuid
print(prepared_protein_uuid)
```

## Settings

- `add_missing_method` (default `"boltz_2"`): use `"boltz_2"` or `"pdbfixer"` to add missing atoms and residues; `None` skips this step.
- `cap_residues` (default `"ace_nme"`): use ACE/NME caps or `"terminal_templates"`; `None` disables capping. ACE/NME requires an add-missing method and cannot be combined with `protonation_method="protonate_utils"`.
- `protonation_method` (default `"openmm"`): use `"openmm"`, `"protonate_utils"`, or `"propka_3"` to add hydrogens.
- `pH` (default `7.4`): pH used to assign protonation states.
- `retain_protonation` (default `True`): retain existing protonation states where possible.
- `retain_non_polymer` (default `{"NA": None, "CL": None, "MG": None}`): non-polymer residues to preserve, keyed by residue name or 0-based residue index. Map ligands and other non-ion residues to their SMILES for parameterization. Known ions and waters may map to `None`. Pass `None` to remove all non-polymer residues.

Example retaining a ligand and a structural water:

```python
wf = rowan.submit_protein_preparation_workflow(
    protein=protein.uuid,
    retain_non_polymer={"LIG": "CC(=O)O", "HOH": None},
)
```

Stored-protein workflow submissions accept any `rowan.Protein` or `rowan.ProteinUUID`. `ProteinUUID` is a semantic string alias: the API verifies that it identifies an accessible protein record. Prefer passing `prepared_protein_uuid` directly when chaining from protein preparation because it avoids fetching and transferring structure data that the local code does not use. Call `get_prepared_protein()` when you need to inspect, modify, or download the prepared structure; it makes one API call on first access, returns a fully loaded `rowan.Protein`, and caches it.

## Fast protein prep

Use `protein.prepare()` when turnaround matters more than the full protein preparation feature set. Protein prep typically finishes in about a minute or less, while full protein preparation can take around ten minutes depending on the structure and settings. It runs PDBFixer and OpenMM, modifies the existing protein record, and blocks until the operation finishes. It can repair missing residues and atoms, remove heterogens, add hydrogens at a selected pH, and optimize hydrogen positions.

```python
protein = rowan.create_protein_from_pdb_id("1CRN")
protein.prepare(add_hydrogen_ph=7.4)
```

Use `submit_protein_preparation_workflow()` when you need Boltz-2 missing-structure modeling, terminal capping, alternative protonation methods, explicit retained non-polymer handling, or an immutable input with a separate prepared output. Use `protein.prepare()` for the faster in-place path when its smaller feature set is sufficient.

## validate_protein_forcefield()

Server-side check that the protein can be parameterized by the MD forcefield. Call before any MD workflow (protein MD, pose-analysis MD, RBFE perturbation) to catch parameterization issues early.

```python
protein.validate_protein_forcefield(exclude_residues=None)
```

Ligand residues (`LIG`) are always excluded automatically. Pass `exclude_residues` for other residues to skip: a residue name excludes the first residue with that name, so further copies are still validated, while an integer is a 0-based index into the protein's sorted non-polymer records and excludes that record without naming it. Both spellings match `rowan.Binder`'s `small_molecules` keys, so `exclude_residues=list(binder.small_molecules)` skips exactly the residues an MD workflow parameterizes from SMILES.

```python
protein.validate_protein_forcefield(exclude_residues=["STI", 2])
```

## Troubleshooting

The preferred preparation workflow strips and reassigns hydrogens automatically. If forcefield validation still fails, inspect the validation error for an unsupported retained cofactor, metal, glycan, or clashing residue. Correct `retain_non_polymer` or the input structure, submit a new preparation workflow, and validate the newly returned protein.

## Pattern by use case

- **Docking, batch docking, analogue docking**: submit protein preparation with defaults, then pass the returned `prepared_protein_uuid`.
- **Protein MD, pose-analysis MD, RBFE perturbation**: retain any bound ligand by mapping its residue name or index to its SMILES in `retain_non_polymer`, then pass `prepared_protein_uuid`. These submitters validate forcefield compatibility by default without downloading the structure.
