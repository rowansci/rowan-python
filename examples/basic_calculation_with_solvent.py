from stjames import Method, Molecule

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")


def compute_energy_with_solvent_correction(molecule: Molecule, method: Method, name: str) -> float:
    opt_workflow = rowan.submit_workflow(
        initial_molecule=molecule,
        workflow_type="basic_calculation",
        name=f"{name} {method} optimization",
        folder_uuid=folder,
        workflow_data={
            "settings": {"method": method, "tasks": ["optimize"]},
        },
    )

    print(f"View workflow privately at: https://labs.rowansci.com/calculation/{opt_workflow.uuid}")
    opt_result = opt_workflow.result()

    sp_workflow = rowan.submit_workflow(
        initial_molecule=opt_result.molecule,
        workflow_type="basic_calculation",
        name=f"{name} {method} single point",
        folder_uuid=folder,
        workflow_data={
            "settings": {
                "method": method,
                "tasks": ["energy"],
                "solvent_settings": {"solvent": "water", "model": "cpcmx"},
            },
        },
    )

    print(f"View workflow privately at: https://labs.rowansci.com/calculation/{sp_workflow.uuid}")
    sp_result = sp_workflow.result()

    return sp_result.energy


E1 = compute_energy_with_solvent_correction(
    Molecule.from_smiles("CC(=O)Oc1ccccc1C(=O)O"), Method.OMOL25_CONSERVING_S, "aspirin"
)
print(E1)
