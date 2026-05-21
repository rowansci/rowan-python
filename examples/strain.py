import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

# Defaults applied (override via `conf_gen_settings=` / `multistage_opt_settings=`):
#   conf_gen: OpenConf, max_confs=200
#   opt: GFN2-xTB / ALPB(water)
#   singlepoint: g-xTB / CPCMx(water)
workflow = rowan.submit_strain_workflow(
    initial_molecule=rowan.Molecule.from_smiles("CCCCCC"),
    name="Hexane Strain",
    folder=folder,
)

print(f"View workflow privately at: https://labs.rowansci.com/strain/{workflow.uuid}")
result = workflow.result()
print(result)
# e.g. <StrainResult strain=2.34 kcal/mol>
