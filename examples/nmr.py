from stjames import Molecule

import rowan

# Set ROWAN_API_KEY environment variable to your API key or set rowan.api_key directly
# rowan.api_key = "rowan-sk..."

workflow = rowan.submit_nmr_workflow(
    Molecule.from_smiles("O[C@H]1[C@H](C(C)C)CC[C@@H](C)C1"),
    name="menthol NMR",
)
print(f"View nmr workflow privately at: https://labs.rowansci.com/workflow/{workflow.uuid}")
workflow.wait_for_result().fetch_latest(in_place=True)

# print hydrogen peaks
for peak in workflow.data["predicted_peaks"]["1"]:
    print(peak)
