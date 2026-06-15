# Relative binding free energy perturbation

## Input

A completed perturbation graph and a prepared protein.

- `graph_result`: a `RelativeBindingFreeEnergyGraphResult` from the RBFE graph workflow, which turns a congeneric ligand series into the graph of ligand pairs FEP transforms between.
- `protein`: a `rowan.Protein` or its UUID. Upload your own PDB with `rowan.upload_protein(name, path)`, or get one from the PDB with `rowan.create_protein_from_pdb_id(name, pdb_code, project_uuid=...)`. Call `protein.prepare()` first to fix nonstandard residues, add missing atoms, and add hydrogens.

This workflow runs FEP simulations along the graph edges to predict relative binding free energies across the ligand series.

## Practical notes

- Accuracy depends heavily on the quality of the input poses. Start from a crystallographic protein–ligand complex when one is available; docked or co-folded poses also work, but only if they are physically valid and free of steric clashes.
- Simulate each ligand in the protonation state and tautomer it actually adopts in solution. Modeling the wrong state can give completely wrong predictions, so check ionizable ligands with a pKa/tautomer tool (e.g. the macro-pKa workflow) before building the graph.
- Expect roughly 10–20 minutes of compute per graph edge with the `recommended` profile. `rigorous` is much slower; `fast` is quicker, though less dramatically so.
- Rowan's FEP engine (TMD) gets its speed largely from *local resampling* — concentrating MD effort on the region around the perturbed atoms rather than the whole system. This is the main lever behind the profiles: `recommended` uses it, while `rigorous` turns it off for an apples-to-apples comparison with conventional FEP.

## Example

```python
from pathlib import Path
import rowan

folder = rowan.get_folder("examples")
data_dir = Path("examples/data")

protein = rowan.upload_protein("TYK2", data_dir / "tyk2_structure.pdb")
protein.prepare()

# Step 1: build the perturbation graph from ligands with 3D coordinates.
ligands = rowan.load_named_ligands(data_dir / "tyk2_ligands.sdf")
graph_wf = rowan.submit_relative_binding_free_energy_graph_workflow(
    ligands=ligands,
    folder=folder,
)
graph_result = graph_wf.result()

# Step 2: run the FEP simulation.
wf = rowan.submit_relative_binding_free_energy_perturbation_workflow(
    graph_result=graph_result,
    protein=protein,
    tmd_settings="recommended",
    folder=folder,
)

result = wf.result()
print(result)   # e.g. <RelativeBindingFreeEnergyPerturbationResult ligands=16>

# Per-ligand binding free energy (dG), relative to the series reference.
for name, res in result.ligand_dg_results.items():
    print(f"  {name}: dG = {res.dg:.2f} +/- {res.dg_err:.2f} kcal/mol")
```

## Settings

The `tmd_settings` profile is the intended interface: it sets sensible defaults for all the lower-level sampling parameters. Pick one and, in most cases, change nothing else.

- `tmd_settings` (default `recommended`): the accuracy/speed profile.
  - `fast`: fewest windows and lowest required overlap, for maximum throughput.
  - `recommended`: the default; uses local resampling to dramatically accelerate the calculation with minimal accuracy loss.
  - `rigorous`: disables local resampling, for careful benchmarking and apples-to-apples comparison with other FEP software.
- `charge_method` (default profile-dependent): partial-charge method. `nagl` is much faster and is recommended whenever any ligand has more than ~50–70 atoms; `amber_am1bcc` is the choice for rigorous calculations.
- `forcefield` (default `off_sage_2_0_0`): simulation force field, `off_sage_2_0_0` or `off_sage_2_2_1`.
- `legs` (default all): which thermodynamic-cycle legs to run, a subset of `vacuum`, `solvent`, `complex`.
- `save_trajectories` (default `False`): save DCD trajectories.
- `trajectory_save_interval` (default `1000`): when saving trajectories, save a frame every N production frames.
- `validate_forcefield` (default `True`): validate protein forcefield compatibility before submitting.

### Advanced sampling parameters

These are FEP internals that the `tmd_settings` profile already sets; override one only with a specific reason. Run `help(rowan.submit_relative_binding_free_energy_perturbation_workflow)` for the rest.

- `n_eq_steps`: equilibration steps before production sampling.
- `n_frames`: total production frames to sample. Each frame is 400 MD steps (`steps_per_frame`), about 1 ps.
- `n_windows`: starting number of lambda windows; the scheduler adaptively reschedules them to cut cost.
- `min_overlap` and `target_overlap`: minimum-required and target overlap between adjacent lambda windows.
- `local_md_steps`: local MD steps per frame, out of the 400; set to `0` to disable local resampling (see Practical notes).
- `local_md_radius` (default `1.2`): radius in nanometers of the local-MD sphere.
- `rest_max_temperature_scale` (default `1.0`): maximum temperature scale for REST (Replica Exchange with Solute Tempering) enhanced sampling; `1.0` disables REST.

## Result fields

- `ligand_dg_results`: per-ligand binding free energy results (`res.dg`, `res.dg_err`), as shown above.
- `edges`: graph edges with per-edge FEP results.
- `diagnostics`: aggregate QC metrics from the FEP simulation.
- `ligands`: the input ligand molecules, keyed by identifier.
