"""
Run a periodic DFT energy calculation on bulk aluminium using Quantum ESPRESSO.

Key settings for PBC calculations:
- ``pw_cutoff``: plane-wave kinetic-energy cutoff in Hartree (higher = more accurate/slower)
- ``kpoints``: Monkhorst–Pack k-point grid (denser = more accurate/slower)
- ``smearing``: occupation smearing type — recommended for metals to aid SCF convergence
- ``degauss``: smearing width in Hartree (typical range: 0.005–0.02)

Periodic molecules are constructed from atomic positions + lattice vectors via stjames.
See documentation at: https://docs.rowansci.com/science/workflows/basic-calculation
"""

import stjames
import stjames.molecule
from stjames.periodic_cell import PeriodicCell

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples/periodic_dft")

# Build bulk Al FCC primitive cell.
# Lattice vectors in Angstrom; Al has 13 electrons so multiplicity=2.
cell = PeriodicCell(
    lattice_vectors=(
        (0.0,    2.0230, 2.0230),
        (2.0230, 0.0,    2.0230),
        (2.0230, 2.0230, 0.0),
    )
)
stj_mol = stjames.Molecule(
    atoms=[stjames.molecule.Atom(atomic_number=13, position=(0.0, 0.0, 0.0))],
    charge=0,
    multiplicity=2,
    cell=cell,
)
al_fcc = rowan.Molecule.from_stjames(stj_mol)

# Marzari–Vanderbilt cold smearing is recommended for metals.
pbc_settings = rowan.PBCDFTSettings(
    pw_cutoff=7.5,           # Hartree; SSSP efficiency recommends ~7–9 Ha for Al
    kpoints=(4, 4, 4),       # Monkhorst–Pack grid; increase for production runs
    smearing=rowan.PBCDFTSmearing.MV,
    degauss=0.01,            # Hartree smearing width
)

workflow = rowan.submit_basic_calculation_workflow(
    initial_molecule=al_fcc,
    tasks=["energy"],
    method="PBE",
    basis_set="SSSP_efficiency",
    pbc_dft_settings=pbc_settings,
    name="Al FCC bulk energy",
    folder=folder,
)

print(f"View workflow privately at: https://labs.rowansci.com/calculation/{workflow.uuid}")
result = workflow.result()
print(result)
# e.g. <BasicCalculationResult energy=-19.725 H>
print(f"Al FCC energy: {result.energy:.6f} Hartree")
