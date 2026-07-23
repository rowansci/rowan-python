# Electronic properties

## Input

A 3D structure (`StructureInput`): a `rowan.Molecule`, `stjames.Molecule`, or RDKit `Mol` carrying coordinates. This workflow does **not** optimize the geometry, so provide an already-optimized structure, for example one loaded with `rowan.Molecule.from_xyz_file(path)`. A raw SMILES embedding (`rowan.Molecule.from_smiles(...)`) is not optimized and will give unreliable properties. If you do not have an optimized structure, run a basic calculation `optimize` first and feed its `result.molecule` here.

## Example

```python
import rowan

folder = rowan.get_folder("examples")

mol = rowan.Molecule.from_xyz_file("formaldehyde_opt.xyz")  # an already-optimized geometry

wf = rowan.submit_electronic_properties_workflow(
    initial_molecule=mol,
    method="b97_3c",
    folder=folder,
)

result = wf.result()
print(result)  # orbitals, density/ESP cubes, charges, bond orders, multipole moments
```

## Settings

- `method` (default `b97_3c`): the level of theory for the wavefunction. It must be a wavefunction method (DFT or Hartree-Fock); neural network potentials and force fields cannot produce orbitals, density, or ESP, so they do not apply. `b97_3c` is a composite method that bundles its own basis set and corrections (a solid lightweight default), so leave `basis_set` as none when using it.
- `basis_set` (default none): set this only for a plain DFT functional that needs one (for example `def2-TZVP` with `pbe0`). Composite methods like `b97_3c` already include a basis set, so leave it none with them.
- `compute_density_cube` and `compute_electrostatic_potential_cube` (default `True`): write volumetric cube grids of the electron density and the electrostatic potential, for visualization (such as mapping ESP onto a surface). These are large; set them `False` if you only need the scalar properties.
- `compute_num_occupied_orbitals` and `compute_num_virtual_orbitals` (default `1`): how many frontier orbitals to export as cubes, counting down from the HOMO (HOMO, HOMO-1, ...) and up from the LUMO (LUMO, LUMO+1, ...). The default exports just the HOMO and LUMO.

Charges (Mulliken, Löwdin), bond orders (Wiberg, Mayer), and multipole moments (dipole, quadrupole) are always computed; only the cubes and orbital exports are gated by the settings above. Both charge schemes and both bond-order schemes are returned: prefer Löwdin charges over Mulliken (better handling of basis-set overlap) and Mayer bond orders over Wiberg (especially for open-shell species).
