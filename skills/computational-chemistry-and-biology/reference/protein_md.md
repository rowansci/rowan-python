# Protein MD

## Input

A `rowan.Protein` (or its UUID string). Proteins are not SMILES. Load one and prepare it before submitting:

- From a PDB ID: `rowan.create_protein_from_pdb_id("1CRN", name="crambin")`.
- From a local PDB file: `rowan.upload_protein("my protein", "path/to/file.pdb")`.

Protein MD accepts any stored `rowan.Protein` or protein UUID. Protein preparation is recommended before submission; when chaining from it, passing `prepared_protein_uuid` avoids an unnecessary structure fetch.

Runs a molecular dynamics simulation on the protein.

## Example

```python
import rowan

folder = rowan.get_folder("examples")

protein = rowan.create_protein_from_pdb_id("1CRN", name="crambin")
preparation_workflow = rowan.submit_protein_preparation_workflow(
    protein=protein.uuid,
    folder=folder,
)
prepared_protein_uuid = preparation_workflow.result().prepared_protein_uuid

wf = rowan.submit_protein_md_workflow(
    protein=prepared_protein_uuid,
    num_trajectories=1,  # example uses 1 for speed; default 4
    simulation_time_ns=1,  # example uses 1 for speed; default 10
    folder=folder,
)

result = wf.result()
print(result.trajectory_uuids)  # one UUID per trajectory
result.download_trajectories([0], path=".")  # save DCD trajectory files
```

## Settings

- `num_trajectories` (default `4`): number of independent trajectories (replicates) to run. More replicates improve conformational sampling; reduce for a quick, cheaper run.
- `small_molecule_ff` (default `off_sage_2_3_0`): force field for small molecules in the structure. `mango_1_0_0` generates ligand-specific parameters using machine learning. OpenFF Sage 2.0.0, 2.2.1, and 2.3.0 are also available.
- `protein_ff` (default `rowan.ProteinForceField.FF14SB`): force field for the protein. Pass a `rowan.ProteinForceField` value.
- `water_ff` (default `rowan.WaterForceField.TIP3P`): force field for water. Pass a `rowan.WaterForceField` value compatible with the protein force field.
- `equilibration_time_ns` (default `0.5`): equilibration time per trajectory, in nanoseconds.
- `simulation_time_ns` (default `10`): production simulation time per trajectory, in nanoseconds.
- `temperature` (default `300`): temperature, in kelvin.
- `pressure_atm` (default `1.0`): pressure, in atmospheres.
- `langevin_timescale_ps` (default `1.0`): timescale for the Langevin integrator, in inverse picoseconds.
- `timestep_fs` (default `4`): integration timestep, in femtoseconds.
- `hydrogen_mass` (default `3`): hydrogen mass in atomic mass units, repartitioned to support the 4 fs timestep.
- `constrain_hydrogens` (default `True`): use SHAKE to freeze bonds to hydrogen.
- `nonbonded_cutoff` (default `8.0`): nonbonded cutoff for particle-mesh Ewald, in angstrom.
- `ionic_strength_M` (default `0.0`): ionic strength of the solution, in molar.
- `water_buffer` (default `8.0`): amount of water added around the protein, in angstrom.
- `save_solvent` (default `False`): whether to save solvent atoms in the trajectories.
- `num_solvent_to_save` (default `None`): when `save_solvent=True` and a `binder` is set, keep only the N solvent molecules nearest the binder each frame; `None` keeps all solvent. Ignored when `save_solvent=False` or no `binder`.
- `small_molecules` (default `None`): SMILES keyed by non-polymer residue name or zero-based index. Use this to parameterize one or more small molecules in the protein; a `None` value selects an existing residue template.
- `binder` (default `None`): a `rowan.Binder` specifying the binder within the complex — protein/peptide chains (`chain_ids`), small molecules (`small_molecule_residues`, identified by residue-name string or 0-based non-polymer residue index), or both. Enables per-frame MM/GBSA and binder RMSD analyses (see result fields).
- `protein_restraint_cutoff` (default `None`): distance from the binder past which Cα atoms are harmonically restrained, in angstrom; `None` disables restraints. Useful for keeping the binding site mobile while stabilizing the rest of the protein.
- `protein_restraint_constant` (default `100`): force constant for the Cα backbone restraints, in kcal/mol/Å².
- `analysis_interval_ps` (default `None`): interval at which to compute per-frame SASA and polar SASA, in ps. `None` (the default) disables those analyses.
- `clustering` (default `None`): cluster the trajectory frames. `None` disables it; pass `rowan.KMeansClusteringSettings(num_clusters=10)` or `rowan.GreedyClusteringSettings(cutoff_angstrom=2.0)`.
- `validate_forcefield` (default `True`): validate the protein forcefield before submitting; raises early if the protein cannot be parameterized or has clashing residues. When a `binder` is set, its small molecules are excluded from validation, whether keyed by residue name or by non-polymer index, since they are parameterized from their SMILES. Cofactors, metals, and glycans outside the binder are still validated — set `False` to skip the pre-check entirely.

## Result fields

- `trajectory_uuids`: UUIDs of the trajectory calculations, one per replicate.
- `trajectories`: per-replicate results. Each exposes `protein_rmsd` (per-frame Cα RMSD from frame 0, Å), `rmsf` (per-Cα RMSF from the mean structure, Å), `potential_energy` (per-frame whole-system potential energy, Hartree), and the radius of gyration per frame (`isotropic_radius_of_gyration`). `sasa` and `polar_sasa` are populated when `analysis_interval_ps` is set; `cluster_centroid_indices` / `cluster_indices_by_frame` when `clustering` is set; and `mmgbsa_scores` (per-frame MM/GBSA binding-side interaction energy, kcal/mol) plus `binder_rmsd` when a `binder` is set.
- `minimized_protein_uuid` / `get_minimized_protein()`: the energy-minimized protein.
- `get_mean_structure(replicate=0)` / `download_mean_structure(...)`: retrieve or download the coordinate-averaged structure for one replicate.
- `download_medoid_structure(replicate=0, ...)`: download the actual trajectory frame closest to the replicate's average structure.
- `bonds`: bond list for the simulated system.
- `messages`: messages emitted during the run.
- `download_trajectories(replicates, path=...)`: download DCD trajectory files for the given replicate indices as a `.tar.gz`.
- `get_atom_distances(atom_pairs, replicate=0)`: fetch per-frame interatomic distances (Angstrom) for a list of `(atom_i, atom_j)` index pairs over the trajectory. Returns one list of floats per pair.
