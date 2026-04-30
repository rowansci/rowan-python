"""
Perform a multistage geometry optimization using the Rowan API.

Defaults to a `r2scan_3c//gfn2_xtb` stack. Pass explicit `optimization_settings`
(a list of `Settings`, one per stage) and `singlepoint_settings` to customize.

See documentation at: https://docs.rowansci.com/science/workflows/multistage-optimization
"""

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

workflow = rowan.submit_multistage_optimization_workflow(
    initial_molecule=rowan.Molecule.from_smiles("C1CCC1"),  # cyclobutane
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
