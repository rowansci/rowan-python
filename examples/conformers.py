"""
Calculate the conformers of a molecule using the Rowan API.

Conformer generation defaults to OpenConf; pass `conf_gen_settings=` to choose a
different generator (ETKDG, iMTD-GC, MCMM).

See documentation at: https://docs.rowansci.com/science/workflows/conformers
"""

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

workflow = rowan.submit_conformer_search_workflow(
    initial_molecule=rowan.Molecule.from_smiles("CCOCC"),
    folder=folder,
)

print(f"View workflow privately at: https://labs.rowansci.com/conformer-search/{workflow.uuid}")

result = workflow.result()
print(result)
# e.g. <ConformerSearchResult conformers=12>
