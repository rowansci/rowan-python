# Covalent inhibitor scan

## Input

A `rowan.Protein` (or its UUID string) with the protein and ligand as separate entities — the ligand as a non-polymer residue. The complex can be covalently bonded already (e.g. from a crystal structure) or just a normally docked (non-covalent) pose; the scan runs from `scan_start` to `scan_stop` regardless of which end the input starts near. Plus two atom indices, both 0-based and in PDB record order (protein atoms, then ligand atoms):

- `protein_reactive_atom_index`: the reacting protein atom (e.g. a catalytic Cys `SG`).
- `ligand_reactive_atom_index`: the reacting ligand atom (e.g. the warhead carbon).

The ligand must be classified as non-polymer in the protein's data — `protein.prepare()` can misclassify a covalently-bonded ligand as part of the polymer chain, in which case the workflow fails at compute time with "Complex PDB has no non-polymer atoms to use as the ligand".

## Example

```python
import rowan

folder = rowan.get_folder("examples")

wf = rowan.submit_covalent_inhibitor_scan_workflow(
    protein=protein,  # covalently bonded or normally docked complex
    protein_reactive_atom_index=1571,
    ligand_reactive_atom_index=4492,
    settings=rowan.CovalentInhibitorScanSettings(scan_num=4),
    folder=folder,
)

result = wf.result()
for distance, energy in result.get_energies():
    print(f"{distance:.3f} Å  {energy}")
```

## Settings

`rowan.CovalentInhibitorScanSettings`:

- `truncation_radius` (default `6.0`): radius in angstrom beyond which the protein is truncated around the ligand.
- `unfreeze_radius` (default `0.0`): side chains with any atom within this radius (angstrom) of any ligand atom are relaxed during the scan; `0` freezes all side chains.
- `scan_start` (default `1.8`), `scan_stop` (default `3.5`): bond-distance range for the scan, in angstrom. Scans covalent to non-covalent by default, but the direction is just whichever end the input geometry is closer to. These defaults are tuned for a C–S bond (e.g. Cys `SG` to an acrylamide warhead carbon) — for other nucleophile/warhead pairs (Ser/Thr `O`, Lys `N`, His `N`, etc.), adjust both to bracket that bond's actual covalent and non-covalent distances.
- `scan_num` (default `10`): number of points along the scan.
- `calc_settings` (default `uma_s_omol`/`omol25`): settings for each scan-point optimization.
- `singlepoint_settings` (default `g_xtb`/`alpb(water)`): settings for the per-point implicit-solvent single-point energy.

## Result fields

- `scan_points`: list of `CovalentInhibitorScanPoint`, each with `index`, `distance` (Å), `molecule`, `energy` (Hartree, or `None` if the single-point failed), and `uuid`.
- `get_energies()`: list of `(distance, energy)` tuples, ordered by scan point index (`scan_start` to `scan_stop`).
