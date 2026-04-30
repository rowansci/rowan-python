"""
Calculate spin state energies using the Rowan API.

This workflow predicts the lowest energy spin state by running multistage
optimizations at different spin multiplicities. The same
`MultiStageOptSettings` is applied to each multiplicity.

See documentation at: https://docs.rowansci.com/science/workflows/spin-states
"""

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

# r2scan_3c//gfn2_xtb stack: GFN2-xTB opt, R²SCAN-3c singlepoint.
mso_settings = rowan.MultiStageOptSettings(
    optimization_settings=[
        rowan.Settings(method=rowan.Method.GFN2_XTB, tasks=[rowan.Task.OPTIMIZE]),
    ],
    singlepoint_settings=rowan.Settings(method=rowan.Method.R2SCAN3C, tasks=[rowan.Task.ENERGY]),
)

# Methylene (CH2) - classic spin states example, triplet is ground state
workflow = rowan.submit_spin_states_workflow(
    initial_molecule=rowan.Molecule.from_smiles("[CH2]"),
    states=[1, 3],  # singlet vs triplet
    multistage_opt_settings=mso_settings,
    name="Methylene spin states",
    folder=folder,
)

print(f"View workflow privately at: https://labs.rowansci.com/spin-states/{workflow.uuid}")

result = workflow.result()
print(result)
# e.g. <SpinStatesResult states=2 ground=(mult=3, E=-39.123456 Ha)>
