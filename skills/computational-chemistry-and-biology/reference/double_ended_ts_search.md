# Double-ended TS search

## Input

Two 3D structures, a `reactant` and a `product`, each a `StructureInput` (a `rowan.Molecule`, `stjames.Molecule`, or RDKit `Mol` carrying coordinates). Get them any way: load coordinates with `rowan.Molecule.from_xyz(...)` or `rowan.Molecule.from_xyz_file(path)`, or embed from a SMILES with `rowan.Molecule.from_smiles(...)`. They must share the same atom ordering so the two endpoints can be interpolated. If a reaction involves multiple reactants or products, combine them into a single 3D structure for that endpoint. If you know the mechanism, arranging the molecules in the reaction channel can speed up the search, but it is not required. The workflow finds the transition state connecting them using a string method (freezing string method, FSM, by default; growing string method, GSM, via `freeze=False`) or a nudged elastic band (NEB).

## Example

```python
import rowan

folder = rowan.get_folder("examples")

reactant = rowan.Molecule.from_xyz("H 0 0 -1.1\nC 0 0 0\nN 0 0 1.2")  # HCN
product = rowan.Molecule.from_xyz("H 0 0 2.3\nC 0 0 0\nN 0 0 1.2")    # CNH

wf = rowan.submit_double_ended_ts_search_workflow(
    reactant=reactant,
    product=product,
    calculation_settings=rowan.Settings(method=rowan.Method.GFN2_XTB),  # level of theory for the search
    folder=folder,
)

result = wf.result()
print(result)   # <DoubleEndedTSSearchResult ts_uuid=... fwd=... bwd=...>
print([p.distance for p in result.forward_path])
print([p.distance for p in result.backward_path])
```

`result.ts_molecule` is the located transition state (with `result.ts_energy`), and `result.get_path_energies(relative=True)` gives the reaction-path energy profile. A double-ended search only produces a TS guess, so verify it actually connects your reactant and product by running an IRC from `result.ts_molecule` (see IRC); this requires `optimize_ts=True` so the TS is a real saddle point.

## Settings

- `calculation_settings` (default none, resolves to omol25_conserving_s): a `rowan.Settings` specifying the level of theory (method, basis set, and so on) used during the search and optimizations. Set an implicit solvent here through its `solvent_settings` (default none, gas phase); see the basic calculation reference for `solvent_settings`.
- `search_settings` (default none, resolves to `rowan.StringMethodSettings()`): either a `rowan.StringMethodSettings` (`freeze=True` -> FSM; `freeze=False` -> GSM) or a `rowan.NEBSettings` (nudged elastic band). Both take `interpolation_method`, a `rowan.Interpolation` (`GEODESIC` by default, or `IDPP`, `CARTESIAN`, `LINEAR_SYNCHRONOUS_TRANSIT`, `REDUNDANT_INTERNAL_COORDINATES`).
- `optimize_inputs` (default `True`): pre-optimize the reactant and product before the search. Leave it on for most cases; turn it off when you have manually arranged the endpoints in the reaction channel, since optimization would relax them out of that arrangement.
- `optimize_ts` (default `True`): optimize the located TS guess to a true transition state.

Double-ended TS searches fail fairly often, which is normal for the method. If a search fails, retry with different endpoints or settings, or email support@rowansci.com.
