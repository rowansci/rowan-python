"""
Perform a multistage geometry optimization using the Rowan API.

Each entry of `optimization_settings` runs in order; `singlepoint_settings`
runs last on the final geometry.

See documentation at: https://docs.rowansci.com/science/workflows/multistage-optimization
"""

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

# r2scan_3c//gfn2_xtb//gfn_ff stack: GFN-FF pre-opt, GFN2-xTB opt, R²SCAN-3c singlepoint.
optimization_settings = [
    rowan.Settings(method=rowan.Method.GFN_FF, tasks=[rowan.Task.OPTIMIZE]),
    rowan.Settings(method=rowan.Method.GFN2_XTB, tasks=[rowan.Task.OPTIMIZE]),
]
singlepoint_settings = rowan.Settings(method=rowan.Method.R2SCAN3C, tasks=[rowan.Task.ENERGY])

workflow = rowan.submit_multistage_optimization_workflow(
    initial_molecule=rowan.Molecule.from_smiles("C1CCC1"),  # cyclobutane
    optimization_settings=optimization_settings,
    singlepoint_settings=singlepoint_settings,
    name="Multistage optimization cyclobutane",
    folder=folder,
)

print(f"View workflow privately at: https://labs.rowansci.com/multistage-opt/{workflow.uuid}")
print(f"Workflow UUID: {workflow.uuid}")

# Stream stage-by-stage progress; final iteration is the complete result.
for result in workflow.stream_result(poll_interval=10):
    print(f"  {len(result.calculation_uuids)} stages done, energy={result.energy}")

print(result)
# e.g. <MultiStageOptResult energy=-76.234567 H>

# Or pick up a completed workflow later by UUID:
# result = rowan.retrieve_workflow(workflow.uuid).result(wait=False)
