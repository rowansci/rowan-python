from stjames import Molecule

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

workflow = rowan.submit_nmr_workflow(
    Molecule.from_smiles("O[C@H]1[C@H](C(C)C)CC[C@@H](C)C1"),
    name="menthol NMR",
    folder=folder,
)
print(f"View nmr workflow privately at: https://labs.rowansci.com/nmr/{workflow.uuid}")
result = workflow.result()

# print hydrogen peaks (atomic number 1 = hydrogen)
for peak in result.predicted_peaks[1]:
    print(peak)
