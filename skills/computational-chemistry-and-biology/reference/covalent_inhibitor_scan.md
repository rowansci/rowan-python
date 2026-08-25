# Covalent inhibitor scan

> **Beta:** This workflow is under active development. Its interface and behavior may change.

## Input

A `rowan.Protein` (or its UUID string) with the protein and ligand as separate entities — the ligand as a non-polymer residue. The complex can be covalently bonded already (e.g. from a crystal structure) or just a normally docked (non-covalent) pose. Plus two atom indices, both 0-based and in PDB record order (protein atoms, then ligand atoms):

- `protein_reactive_atom_index`: the reacting protein atom (e.g. a catalytic Cys `SG`).
- `ligand_reactive_atom_index`: the reacting ligand atom (e.g. the warhead carbon).
- `reactant_smiles`: SMILES of the neutral reactant ligand, used to rebuild the ligand as a separate non-covalent molecule at the MM level.

The ligand must be classified as non-polymer in the protein's data. When preparing the complex, pass the ligand residue name or index and its SMILES through `retain_non_polymer`. Verify the ligand remains a non-polymer in the prepared output; a covalently bonded ligand classified as part of the polymer chain causes the workflow to fail at compute time with "Complex PDB has no non-polymer atoms to use as the ligand". Measure both reactive atom indices from the final prepared structure, because preparation can change atom ordering and residue numbering. Use `protein.get_atom_index(...)` to convert chain, residue, and atom names to the required indices.

If you don't already have a covalently-bonded complex (e.g. from a crystal structure), covalent docking (see `docking.md`, "gnina (noncovalent and covalent docking)") is a good way to generate one: it forms the bond between the specified reactive atoms and returns a complex whose `protein_reactive_atom_index`/`ligand_reactive_atom_index` (there, `covalent_protein_atom_index`/`covalent_ligand_atom_index`) can feed directly into this scan.

The workflow seeds windows along the reactive bond distance with a steered pull, samples each via ML/MM umbrella sampling, and combines them into a free energy profile (not a series of independent single-point energies).

## Example

```python
import rowan

folder = rowan.get_folder("examples")

# Prepared 4YHF: preparation renumbers Cys481 to A.101 and retains ligand 4C9 as A.701.
protein_reactive_atom_index = protein.get_atom_index("A", 101, "SG")
ligand_reactive_atom_index = protein.get_atom_index("A", 701, "C1", entity_type="non_polymer")

wf = rowan.submit_covalent_inhibitor_scan_workflow(
    protein=protein.uuid,  # covalently bonded or normally docked complex
    protein_reactive_atom_index=protein_reactive_atom_index,
    ligand_reactive_atom_index=ligand_reactive_atom_index,
    reactant_smiles=ligand_smiles,
    folder=folder,
)

result = wf.result()
for distance, free_energy in result.get_energies():
    print(f"{distance:.3f} Å  {free_energy}")
print(f"barrier: {result.barrier} kcal/mol")
```

## Settings

`rowan.UmbrellaSamplingScanSettings`:

- `window_centers` (default spans ~1.8–4.0 Å): bias centers along the reactive bond, in angstrom. Tuned for a C–S bond (e.g. Cys `SG` to an acrylamide warhead carbon) — for other nucleophile/warhead pairs (Ser/Thr `O`, Lys `N`, His `N`, etc.), pass a custom list bracketing that bond's actual covalent and non-covalent distances.
- `force_constant` (default `30.0`): umbrella harmonic bias force constant, in kcal/mol/Å².
- `protein_restraint_cutoff` (default `10.0`): distance from the ligand beyond which backbone Cα atoms are harmonically restrained, in Å; `None` disables restraints.
- `protein_restraint_constant` (default `100.0`): Cα restraint force constant, in kcal/mol/Å².
- `calc_settings` (default `orb_v3_conservative_omol`/orb engine): settings for the ML subengine on the model region.
- Also inherits the general protein-MD settings (`equilibration_time_ns`, `simulation_time_ns`, `timestep_fs`, `hydrogen_mass`, `constrain_hydrogens`, forcefields, etc.) — see `protein_md.md`.

## Result fields

- `points`: list of `CovalentInhibitorScanPoint`, each with `index`, `distance` (Å), `force_constant` (kcal/mol/Å²), `free_energy` (kcal/mol, or `None` if unavailable), `mean_distance` (Å), `n_samples`, and `molecule`.
- `convergence`: `UmbrellaSamplingConvergence` diagnostics (`overlap_matrix`, `round_trips`, `worst_pair_acceptance`), or `None`.
- `barrier`: activation free energy for addition, in kcal/mol, or `None` if the profile has no interior maximum.
- `get_energies()`: list of `(distance, free_energy)` tuples, ordered by window index.
