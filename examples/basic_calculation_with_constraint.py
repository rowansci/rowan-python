import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

molecule = rowan.Molecule.from_smiles("CCCC")

# Constrain the C-C-C-C dihedral angle to 0 degrees during optimization.
# Constraint types: "bond", "angle", "dihedral", "freeze_atoms"
# Atoms are 1-indexed.
workflow = rowan.submit_basic_calculation_workflow(
    initial_molecule=molecule,
    method="gfn2_xtb",
    tasks=["optimize"],
    name="Constrained Butane",
    folder=folder,
    opt_settings=rowan.OptimizationSettings(
        constraints=[
            rowan.Constraint(atoms=[4, 3, 2, 1], constraint_type="dihedral", value=0),
        ]
    ),
)

print(f"View workflow privately at: https://labs.rowansci.com/calculation/{workflow.uuid}")
result = workflow.result()
print(result)
# e.g. <BasicCalculationResult energy=-13.657247 H>
