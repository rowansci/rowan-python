"""Band structure of bulk silicon using periodic DFT (Quantum ESPRESSO)."""

import rowan

# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

a = 5.431  # Å
si = rowan.Molecule.from_atoms(
    atoms=[
        rowan.Atom(atomic_number=14, position=(0.0, 0.0, 0.0)),
        rowan.Atom(atomic_number=14, position=(a / 4, a / 4, a / 4)),
    ],
    charge=0,
    multiplicity=1,
    cell=rowan.PeriodicCell(
        lattice_vectors=((0.0, a / 2, a / 2), (a / 2, 0.0, a / 2), (a / 2, a / 2, 0.0))
    ),
)

workflow = rowan.submit_basic_calculation_workflow(
    initial_molecule=si,
    tasks=["band_structure"],
    method="PBE",
    basis_set="SSSP_PBE_efficiency",
    pbc_dft_settings=rowan.PBCDFTSettings(kpoints=(2, 2, 2)),
    name="Si band structure",
    folder=folder,
)
print(f"https://labs.rowansci.com/calculation/{workflow.uuid}")
result = workflow.result()

print(f"Symmetry:  space group {result.symmetry}")
print(f"XRD peaks: {len(result.xrd_peaks)} reflections")
print(f"Band gap:  {result.band_gap:.4f} Ha")  # ~0.020 Ha (PBE underestimates Si's 1.12 eV gap)
print(f"VBM:       {result.band_structure.valence_band_maximum:.4f} Ha")
print(f"CBM:       {result.band_structure.conduction_band_minimum:.4f} Ha")
print(f"DOS:       {len(result.density_of_states)} points")
