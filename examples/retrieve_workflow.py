import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."

# Retrieve a previously submitted workflow by UUID and get its result.
# The UUID can be found in three ways:
#   1. Printed when you submit: print(workflow.uuid)
#   2. The last segment of the workflow URL in the Rowan UI:
#      https://labs.rowansci.com/.../<uuid>
#   3. Listed via rowan.list_workflows()
workflow = rowan.retrieve_workflow("619fedc2-d8e9-4bbc-8780-170108d1f674")

result = workflow.result()
print(result)
# <RelativeBindingFreeEnergyPerturbationResult ligands=16>

print(result.diagnostics)
# RelativeBindingFreeEnergyDiagnostics(
#     cycle_closure_rms=None, windows_completed=52, windows_failed=None
# )

for name, res in result.ligand_dg_results.items():
    print(f"  {name}: ddG = {res.dg:.2f} +/- {res.dg_err:.2f} kcal/mol")
# ejm_31: ddG = 1.79 +/- 0.09 kcal/mol
# ejm_42: ddG = 1.27 +/- 0.08 kcal/mol
# ...

# You can also list recent workflows:
workflows = rowan.list_workflows()
for wf in workflows:
    print(f"  {wf.uuid}  {wf.name}")
