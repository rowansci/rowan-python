# Pose-analysis MD

Useful for analyzing pose stability and specific protein–ligand interactions, but not a replacement for FEP. Effectiveness varies significantly by project; experimental data early and often is advised.

## Input

A docked protein-ligand complex (a holo structure) and the ligand SMILES. The protein must contain the bound ligand as a residue, since the workflow runs MD on the complex and measures the ligand's RMSD against its starting pose; a plain apo protein has no pose to analyze.

- `protein`: the holo complex on which MD runs, as a `rowan.Protein` or its UUID. Upload a complex with `rowan.upload_protein(name, path)`, or obtain one from a co-folding prediction. Either way, call `protein.prepare(remove_heterogens=False)` first to add hydrogens while keeping the bound ligand.
- `initial_smiles`: the SMILES of the ligand bound in the complex, the same molecule. It supplies the chemical identity used to parameterize the ligand for MD (bond orders, force-field atom types); the 3D pose comes from the complex structure.

The MD runs in explicit solvent, and backbone atoms outside the binding pocket are restrained to preserve the protein's fold while the pocket and ligand stay flexible.

## Example

```python
import rowan

folder = rowan.get_folder("examples")

ligand = "CCC(C)(C)NC1=NCC2(CCC(=O)C2C)N1"

# load a docked protein-ligand complex; remove_heterogens=False keeps the bound ligand
protein = rowan.upload_protein("complex", "complex.pdb")
protein.prepare(remove_heterogens=False)

wf = rowan.submit_pose_analysis_md_workflow(
    protein=protein,
    initial_smiles=ligand,
    num_trajectories=1,  # example uses 1 for speed; default 4
    simulation_time_ns=1,  # example uses 1 for speed; default 10
    folder=folder,
)

result = wf.result()
print(result.trajectories[0].ligand_rmsd)  # ligand RMSD by frame
```

## Settings

- `num_trajectories` (default `4`): number of independent replicas to run. More replicas give a more robust read on pose stability; reduce it for a quick, cheaper check.
- `equilibration_time_ns` (default `1`): equilibration time per trajectory, in ns.
- `simulation_time_ns` (default `10`): production simulation time per trajectory, in ns.
- `temperature` (default `300`): temperature, in K.
- `pressure_atm` (default `1.0`): pressure, in atm.
- `langevin_timescale_ps` (default `1.0`): Langevin integrator timescale, in ps^-1.
- `timestep_fs` (default `2`): integration timestep, in femtoseconds.
- `constrain_hydrogens` (default `True`): use SHAKE to freeze bonds to hydrogen.
- `nonbonded_cutoff` (default `8.0`): nonbonded cutoff for particle-mesh Ewald, in angstroms.
- `ionic_strength_M` (default `0.0`): ionic strength of the solution, in molar.
- `water_buffer` (default `10.0`): water padding added around the protein, in angstroms.
- `ligand_residue_name` (default `LIG`): residue name of the ligand in the structure.
- `protein_restraint_cutoff` (default `7.0`): distance from the ligand past which alpha-carbons are restrained, in angstroms. Residues within this cutoff of the pocket move freely. Set to `None` to apply no restraints.
- `protein_restraint_constant` (default `100`): force constant for backbone restraints, in kcal/mol/angstrom^2.
- `save_solvent` (default `False`): save solvent molecules in the trajectory.
- `num_solvent_to_save` (default `None`): when `save_solvent` is on, save only the N solvent molecules nearest the ligand each frame; `None` saves all of them.
- `analysis_interval_ps` (default `None`): interval at which to compute per-frame SASA and polar SASA, in ps. `None` (the default) disables those analyses.
- `clustering` (default `None`): cluster the trajectory frames. `None` disables it; pass `rowan.KMeansClusteringSettings(num_clusters=10)` or `rowan.GreedyClusteringSettings(cutoff_angstrom=2.0)`.
- `validate_forcefield` (default `True`): not a simulation setting but a client-side pre-check. Before submitting, it runs the protein's `validate_protein_forcefield()` (the same call you can run yourself) and raises early if the protein cannot be parameterized or has clashing residues. Set it `False` if you have already validated the protein.

## Result fields

- `trajectories`: per-replicate results. Each exposes the ligand RMSD by frame (`ligand_rmsd`), the persistent protein-ligand contacts with their occupancy (`contacts`), and the radius of gyration per frame (`isotropic_radius_of_gyration`). When `analysis_interval_ps` is set, `sasa` and `polar_sasa` are populated; when `clustering` is set, `cluster_centroid_indices` and `cluster_indices_by_frame` are populated.
- `hydration_sites`: list of `rowan.HydrationSite` objects identified across all trajectories. Each site has a `centroid` (x, y, z tuple in Å), `abs_occ` (absolute occupancy 0–1), `prot_only_abs_occ`, `lig_only_abs_occ`, `bridge_abs_occ`, `num_unique_waters`, and `bridge_residues` (list of `rowan.HydrationBridgeResidue` with `residue_index` and `bridge_occ`).
- `average_rmsds`: average ligand RMSD per trajectory, in angstrom (a low value means the pose held).
- `minimized_protein_uuid` / `get_minimized_protein()`: UUID and fetched `Protein` object of the energy-minimized structure.
- `messages`: messages or warnings emitted during the run.
- `get_atom_distances(atom_pairs, replicate=0)`: fetch per-frame interatomic distances (Angstrom) for a list of `(atom_i, atom_j)` index pairs over the trajectory. Atom indices appear in each trajectory's `contacts` field (`ligand_atom_index`, `protein_atom_index`). Returns one list of floats per pair.
- `download_trajectories(replicates, path=...)`: download DCD trajectory files for the given replicate indices as a `.tar.gz`.
