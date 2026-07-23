# Protein MD

## Input

A `rowan.Protein` (or its UUID string). Proteins are not SMILES. Load one and prepare it before submitting:

- From a PDB ID: `rowan.create_protein_from_pdb_id("crambin", "1CRN", project_uuid=rowan.default_project().uuid)`.
- From a local PDB file: `rowan.upload_protein("my protein", "path/to/file.pdb")`.

Then call `protein.prepare()`, which runs PDBFixer to fix nonstandard residues and add missing atoms and hydrogens, and waits for it to finish.

Runs a molecular dynamics simulation on the protein.

## Example

```python
import rowan

folder = rowan.get_folder("examples")

protein = rowan.create_protein_from_pdb_id(
    "crambin", "1CRN", project_uuid=rowan.default_project().uuid
)
protein.prepare()  # fix residues, add missing atoms/hydrogens; blocks until done

wf = rowan.submit_protein_md_workflow(
    protein=protein,
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
- `equilibration_time_ns` (default `1`): equilibration time per trajectory, in nanoseconds.
- `simulation_time_ns` (default `10`): production simulation time per trajectory, in nanoseconds.
- `temperature` (default `300`): temperature, in kelvin.
- `pressure_atm` (default `1.0`): pressure, in atmospheres.
- `langevin_timescale_ps` (default `1.0`): timescale for the Langevin integrator, in inverse picoseconds.
- `timestep_fs` (default `2`): integration timestep, in femtoseconds.
- `constrain_hydrogens` (default `True`): use SHAKE to freeze bonds to hydrogen.
- `nonbonded_cutoff` (default `8.0`): nonbonded cutoff for particle-mesh Ewald, in angstrom.
- `ionic_strength_M` (default `0.0`): ionic strength of the solution, in molar.
- `water_buffer` (default `10.0`): amount of water added around the protein, in angstrom.
- `save_solvent` (default `False`): whether to save solvent atoms in the trajectories.
- `analysis_interval_ps` (default `None`): interval at which to compute per-frame SASA and polar SASA, in ps. `None` (the default) disables those analyses.
- `clustering` (default `None`): cluster the trajectory frames. `None` disables it; pass `rowan.KMeansClusteringSettings(num_clusters=10)` or `rowan.GreedyClusteringSettings(cutoff_angstrom=2.0)`.
- `validate_forcefield` (default `True`): validate the protein forcefield before running.

## Result fields

- `trajectory_uuids`: UUIDs of the trajectory calculations, one per replicate.
- `trajectories`: per-replicate results. Each exposes the radius of gyration per frame (`isotropic_radius_of_gyration`); `sasa` and `polar_sasa` when `analysis_interval_ps` is set; and `cluster_centroid_indices` / `cluster_indices_by_frame` when `clustering` is set.
- `minimized_protein_uuid` / `get_minimized_protein()`: the energy-minimized protein.
- `bonds`: bond list for the simulated system.
- `messages`: messages emitted during the run.
- `download_trajectories(replicates, path=...)`: download DCD trajectory files for the given replicate indices as a `.tar.gz`.
- `get_atom_distances(atom_pairs, replicate=0)`: fetch per-frame interatomic distances (Angstrom) for a list of `(atom_i, atom_j)` index pairs over the trajectory. Returns one list of floats per pair.
