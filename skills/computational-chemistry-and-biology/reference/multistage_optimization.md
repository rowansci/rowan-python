# Multistage optimization

## Input

A 3D structure (`StructureInput`): a `rowan.Molecule`, `stjames.Molecule`, or RDKit `Mol` carrying coordinates. Get one any way: load coordinates with `rowan.Molecule.from_xyz_file(path)`, reuse a prior result's `.molecule`, or embed from a SMILES with `rowan.Molecule.from_smiles(...)`.

## Example

```python
import rowan

folder = rowan.get_folder("examples")

# omit the stage settings to use the default r2scan_3c//gfn2_xtb stack
wf = rowan.submit_multistage_optimization_workflow(
    initial_molecule=rowan.Molecule.from_smiles("C1CCC1"),
    folder=folder,
)

result = wf.result()
print(result.energy)   # energy in Hartree
```

## Settings

When both `optimization_settings` and `singlepoint_settings` are omitted, the workflow runs a default stack: a GFN2-xTB optimization followed by an r2SCAN-3c single point (`r2scan_3c//gfn2_xtb`). Set them only when you need a different stack. The default suits most work; for publication-quality energetics, use a more accurate stack: a GFN2-xTB optimization, then an r2SCAN-3c optimization, then a ωB97X-3c single point.

- `optimization_settings` (default none): a list of `Settings` stages run in order, each optimizing the previous stage's geometry, e.g. `[rowan.Settings(method=..., tasks=[rowan.Task.OPTIMIZE]), ...]` (Settings are built as in basic calculation). Order them cheap-and-coarse to accurate; a cheap pre-optimization is useful for cleaning up a rough or hand-drawn geometry before the expensive final stage.
- `singlepoint_settings` (default none): a final `Settings` single point run on the optimized geometry.
- `frequencies` (default `False`): also compute vibrational frequencies on the final geometry.
- `transition_state` (default `False`): optimize to a transition state rather than a minimum. Applies to every optimization stage, so a transition-state guess is not relaxed to a minimum along the way. Provide a transition-state guess as the input.
