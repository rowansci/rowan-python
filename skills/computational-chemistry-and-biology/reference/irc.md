# IRC

## Input

A 3D transition-state structure (`StructureInput`: a `rowan.Molecule`, `stjames.Molecule`, or RDKit `Mol` carrying coordinates), typically from a transition-state optimization (the `optimize_ts` task or a double-ended TS search). Load coordinates with `rowan.Molecule.from_xyz_file(path)`, or reuse a prior result's `.molecule` (for example a double-ended TS search's `ts_molecule`). The input must be a transition state, not an arbitrary geometry, so a raw SMILES embedding will not work.

## Example

```python
import rowan

folder = rowan.get_folder("examples")

wf = rowan.submit_irc_workflow(
    initial_molecule=rowan.Molecule.from_xyz(ts_xyz),
    method="omol25_conserving_s",
    folder=folder,
)

result = wf.result()
print(result)  # forward and backward reaction paths

# the two species the TS connects are the ends of each branch:
reactant = result.backward_molecules[-1]
product = result.forward_molecules[-1]
# compare these to your expected reactant/product to confirm the TS connects them
profile = result.get_forward_energies(relative=True)  # kcal/mol along the forward branch
```

## Settings

- `method` (required): level of theory. Ideally use the same level of theory the transition state was optimized at. `basis_set`, `corrections`, and `engine` are also accepted, like a basic calculation.
- `solvent_settings` (default gas phase): implicit solvent (see Solvent below).
- `preopt` (default `True`): re-optimize the input transition state before tracing the path. Use with caution unless the input geometry is already very close to the transition state, since re-optimizing a poor guess can move it to a different saddle point.
- `step_size` (default `0.05`, range `0.001`-`0.5`, in amu^(1/2)·Å) and `max_irc_steps` (default `30`): Initial IRC step size and maximum steps in each direction. The defaults usually suffice; use a smaller `step_size` for a smoother path (especially proton-transfer reactions) and more `max_irc_steps` to capture reactions with large-scale motion. Actual per-step sizes taken are available on the result as `forward_step_sizes`/`backward_step_sizes` once the workflow completes (they are not streamed).
- `optimize_endpoints` (default `False`): optimize the forward and backward endpoint geometries once the IRC completes. When enabled, `result.forward_endpoint_molecule` / `backward_endpoint_molecule` give the optimized product/reactant; use to confirm the endpoints connect to the desired species or to obtain reaction-barrier heights.

## Solvent

Add an implicit solvent by passing `solvent_settings` as a dict (or `rowan.SolventSettings`) with a `solvent` and a `model`, e.g. `{"solvent": "water", "model": "cpcmx"}`. Valid models depend on the engine (an unsupported one raises `ValueError`); check `rowan.ENGINE_SOLVENT_MODELS[engine]`, and `rowan.Solvent` for the available solvents.
