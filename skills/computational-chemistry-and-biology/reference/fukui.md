# Fukui

## Input

A 3D structure (`StructureInput`): a `rowan.Molecule`, `stjames.Molecule`, or RDKit `Mol` carrying coordinates. Get one any way: load coordinates with `rowan.Molecule.from_xyz_file(path)`, reuse a prior result's `.molecule`, or embed from a SMILES with `rowan.Molecule.from_smiles(...)`. The geometry is optimized before the Fukui indices are computed.

## Example

```python
import rowan

folder = rowan.get_folder("examples")

wf = rowan.submit_fukui_workflow(
    initial_molecule=rowan.Molecule.from_smiles("C1=CC=C(C=C1)C(=O)O"),
    optimization_method="gfn2_xtb",
    fukui_method="gfn1_xtb",
    folder=folder,
)

result = wf.result()
print(result)   # per-atom Fukui indices and the global electrophilicity index
```

## Settings

- `optimization_method` (default `gfn2_xtb`): method used to optimize the geometry first.
- `fukui_method` (default `gfn1_xtb`): method used to compute the charges the Fukui indices are derived from.
- `solvent_settings` (default gas phase): implicit solvent (see Solvent below).

The two methods are tuned independently. GFN1-xTB and GFN2-xTB (the defaults) are good for routine, low-cost work; switch `fukui_method` to a DFT method for more accurate charges in difficult cases.

## Solvent

Add an implicit solvent by passing `solvent_settings` as a dict (or `rowan.SolventSettings`) with a `solvent` and a `model`, e.g. `{"solvent": "water", "model": "cpcmx"}`. Valid models depend on the engine (an unsupported one raises `ValueError`); the default GFN methods run on `xtb`, which supports `alpb`, `cpcmx`, and `gbsa`. Check `rowan.ENGINE_SOLVENT_MODELS[engine]`, and `rowan.Solvent` for the available solvents. The geometry optimization always runs in gas phase; the implicit solvent is applied only to the final charge calculation.
