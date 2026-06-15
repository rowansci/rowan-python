import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

# Defaults applied (override via `conf_gen_settings=` / `multistage_opt_settings=`):
#   conf_gen: OpenConf, max_confs=20
#   opt: AIMNet2/wB97M-D3 (gas phase)
#   singlepoint: AIMNet2/wB97M-D3 / CPCMx(water)
workflow = rowan.submit_tautomer_search_workflow(
    initial_molecule=rowan.Molecule.from_smiles("O=C1C=CC=CN1"),
    name="2-Pyridone Tautomers",
    folder=folder,
)

print(f"View workflow privately at: https://labs.rowansci.com/tautomers/{workflow.uuid}")
result = workflow.result()
print(result)
# e.g. <TautomerResult tautomers=3 lowest_energy=-323.456789 H>
