# Basic calculation

## Input

Runs on a single 3D structure (`StructureInput`: a `rowan.Molecule`, `stjames.Molecule`, or RDKit `Mol` carrying coordinates). Get one by:

- generating 3D coordinates from a SMILES string (2D topology that RDKit embeds into 3D): `rowan.Molecule.from_smiles("CCO")` (a bare SMILES string is not accepted)
- reading existing 3D coordinates: `rowan.Molecule.from_xyz(xyz_string)` or `rowan.Molecule.from_xyz_file(path)`
- reusing a structure from a previous result, such as feeding an optimized geometry into a single point: `prev_result.molecule`

## Tasks

Request one or more in a list, for example `["optimize", "frequencies"]`:

- `energy`: single-point energy
- `optimize`: relax to a nearby minimum
- `optimize_ts`: optimize to a transition state. Start from a guess TS structure, not an equilibrium geometry.
- `frequencies`: vibrational frequencies and thermochemistry
- `band_structure`: electronic band structure, DOS, and band gap (periodic systems only)
- `elastic_tensor`: elastic stiffness matrix in GPa (periodic systems only)
- also available: `gradient`, `hessian`, `dipole`, `charge`, `spin_density`, `stress`

Running `frequencies` after an `optimize` or `optimize_ts` is recommended: it confirms whether a local minimum (no imaginary frequencies) or a saddle point (one imaginary frequency) was found, and it computes thermochemical properties. It does increase compute time.

## Presets

A preset bundles a method, engine, basis set, and dispersion corrections. Pass `preset=` together with `tasks`. It is mutually exclusive with `method`, `engine`, `basis_set`, and `corrections` (passing both raises `ValueError`), but you can still pass `mode`, `solvent_settings`, and `opt_settings` alongside it.

- `general_nnp`: omol25_conserving_s on the omol25 engine
- `organic_nnp`: aimnet2_wb97md3 on the aimnet2 engine
- `rapid_semiempirical`: gfn2_xtb on the xtb engine
- `routine_dft`: r2SCAN-D4 / vDZP on the gpu4pyscf engine
- `careful_dft`: wB97M-D3BJ / vDZP on the gpu4pyscf engine

Prefer a preset unless you have a specific reason to set `method`, `engine`, `basis_set`, and `corrections` yourself. Hand-assembled combinations are where invalid or unusual settings come from. If you do set them raw, confirm the combination (see Checking compatibility) or submit and read the `ValueError`.

## Example

```python
import rowan

folder = rowan.get_folder("examples")

wf = rowan.submit_basic_calculation_workflow(
    initial_molecule=rowan.Molecule.from_smiles("CC(=C)C=C"),
    tasks=["optimize"],
    preset="general_nnp",
    folder=folder,
)

result = wf.result()
print(result.energy)  # energy in Hartree
print(result.molecule)  # resulting structure with computed properties
```

## Tuning individual settings

Instead of a preset, set settings directly. Run `help(rowan.submit_basic_calculation_workflow)` for the full signature and defaults. The choices that matter, and the ones agents get wrong:

- `method` (default `omol25_conserving_s`): level of theory, spanning DFT, neural network potentials, semiempirical methods, and force fields. Introspect `rowan.Method` for the full list. A method does not run on every engine, take a basis set, or accept corrections (see Checking compatibility).
- `basis_set`: applies only to DFT methods, e.g. `"def2-TZVP"`. NNP and semiempirical methods ignore it.
- `corrections`: dispersion corrections for DFT, e.g. `["d3bj"]`. NNP and semiempirical methods take none (`gfn2_xtb` with any correction raises `ValueError`).
- `mode` (default `auto`): geometry optimization convergence. Use `careful` for the geometry optimizations you intend to trust and `rapid` for single points or rough/screening passes. `auto` resolves to `rapid`, so a default-mode `optimize` is loosely converged. Note `mode` sets the gradient and energy thresholds and overwrites any you pass in `opt_settings`.
- `solvent_settings` (default gas phase): see Solvent below.
- `opt_settings`: optimization options other than the gradient/energy thresholds (constraints, max steps, transition state). See Optimization settings below.

The named `mode` values set these convergence thresholds:

| Mode | ∆ energy (Hartree) | Max gradient (Hartree/Å) | RMS gradient (Hartree/Å) |
| --- | --- | --- | --- |
| `reckless` | 0.00002 | 0.007 | 0.006 |
| `rapid` | 0.00005 | 0.005 | 0.0035 |
| `careful` | 0.000001 | 0.0009 | 0.0006 |
| `meticulous` | 0.000001 | 0.00003 | 0.00002 |
| `debug` | 0.000001 | 0.000004 | 0.000002 |

## Checking compatibility

Not every method, engine, correction, and basis-set combination is valid. An invalid one raises a `ValueError` at submission naming the supported options, so you can also just submit and read the error. To check first, introspect:

- `rowan.METHOD_ENGINES[method]`: engines that run a method.
- `rowan.ENGINE_METHODS[engine]`: methods an engine supports.
- `rowan.get_supported_corrections(method, engine)`: valid corrections.
- `rowan.ENGINE_SUPPORTS_BASIS_SET[engine]`: whether a basis set applies.
- `rowan.ENGINE_SOLVENT_MODELS[engine]`: valid solvent models.

## Optimization settings

`opt_settings=rowan.OptimizationSettings(...)` controls optimization options that `mode` does not set: constraints, `max_steps`, displacement thresholds, transition-state optimization, and cell relaxation. Run `help(rowan.OptimizationSettings)` for the full list and defaults. Gradient and energy convergence come from `mode`, not here: thresholds you set in `opt_settings` for those are overwritten.

Constraint types are `bond`, `angle`, `dihedral`, and `freeze_atoms`, with atoms 1-indexed:

```python
wf = rowan.submit_basic_calculation_workflow(
    initial_molecule=rowan.Molecule.from_smiles("CCCC"),
    tasks=["optimize"],
    method="gfn2_xtb",
    engine="xtb",
    mode="careful",
    opt_settings=rowan.OptimizationSettings(
        constraints=[
            rowan.Constraint(atoms=[4, 3, 2, 1], constraint_type="dihedral", value=0),
        ],
    ),
    folder=folder,
)
```

## Solvent

Add an implicit solvent by passing `solvent_settings` as a dict (or `rowan.SolventSettings`) with a `solvent` and a `model`:

```python
wf = rowan.submit_basic_calculation_workflow(
    initial_molecule=rowan.Molecule.from_smiles("CC(=O)O"),
    tasks=["optimize"],
    method="gfn2_xtb",
    engine="xtb",
    solvent_settings={"solvent": "water", "model": "cpcmx"},
    folder=folder,
)
```

The valid `model` values depend on the engine, and an unsupported one raises `ValueError`. Check the per-engine list with `rowan.ENGINE_SOLVENT_MODELS[engine]`, and `rowan.Solvent` for the available solvents.

## Periodic systems (PBC)

A periodic system (crystal or surface) needs a unit cell, set either in the loaded structure (an xyz with a `cell:` comment line, or a structure file with lattice vectors) or with `mol.cell = rowan.PeriodicCell(lattice_vectors=...)`. This holds for every engine, including NNPs. Periodic runs work on NNPs (omol25/UMA on `omol25`), GFN (`tblite`), and periodic DFT (`quantum_espresso`); add `opt_settings=OptimizationSettings(optimize_cell=True)` to relax the cell with the atoms.

```python
mol = rowan.Molecule.from_xyz_file("crystal.xyz")  # must contain the unit cell

wf = rowan.submit_basic_calculation_workflow(
    initial_molecule=mol,
    tasks=["optimize"],
    method="uma_s_omat",
    engine="omol25",
    opt_settings=rowan.OptimizationSettings(optimize_cell=True),
    folder=folder,
)
```

For periodic DFT, pass `pbc_dft_settings=rowan.PBCDFTSettings(...)`, which selects the `quantum_espresso` engine (NNPs and GFN do not need it):

```python
wf = rowan.submit_basic_calculation_workflow(
    initial_molecule=mol,
    tasks=["optimize"],
    method="pbe",
    engine="quantum_espresso",
    pbc_dft_settings=rowan.PBCDFTSettings(plane_wave_cutoff=30.0, kpoints=(4, 4, 4)),
    opt_settings=rowan.OptimizationSettings(optimize_cell=True),
    folder=folder,
)
```

`PBCDFTSettings` fields (all optional — auto-derived from pseudopotentials when omitted):
- `plane_wave_cutoff`: kinetic-energy cutoff in Hartree
- `kpoints`: Monkhorst–Pack grid as `(nx, ny, nz)`
- `smearing_type`: `rowan.PBCDFTSmearing.MARZARI_VANDERBILT` (recommended for metals)
- `smearing_width`: smearing width in Hartree (default 0.005)
- `hubbard_u`: DFT+U per element, e.g. `{"Fe": 4.0}`, or `"auto"`

## PBC result properties

Every completed PBC calculation populates these on the result:

- `result.symmetry`: space group number (1–230)
- `result.xrd_peaks`: powder XRD reflections as `(h, k, l, intensity)` tuples

With `tasks=["band_structure"]`:

- `result.band_gap`: band gap in Hartree
- `result.band_structure`: full `BandStructure` object with `eigenvalues`, `kpoint_distances`, `high_symmetry_points`, `valence_band_maximum`, `conduction_band_minimum`
- `result.density_of_states`: total DOS as `(energy, count)` pairs (Fermi level = 0)

With `tasks=["elastic_tensor"]`:

- `result.elastic_tensor`: 6×6 stiffness matrix in GPa (Voigt order: xx yy zz yz xz xy)
