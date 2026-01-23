from stjames import Method, Molecule

import rowan

# Set ROWAN_API_KEY environment variable to your API key or set rowan.api_key directly
# rowan.api_key = "rowan-sk..."


def compute_energy_with_solvent_correction(molecule: Molecule, method: Method, name: str) -> float:
    opt_workflow = rowan.submit_workflow(
        initial_molecule=molecule,
        workflow_type="basic_calculation",
        name=f"{name} {method} optimization",
        workflow_data={
            "settings": {"method": method, "tasks": ["optimize"]},
        },
    )

    print(f"View workflow privately at: https://labs.rowansci.com/workflow/{opt_workflow.uuid}")
    opt_workflow.wait_for_result().fetch_latest(in_place=True)

    calculation_uuid = opt_workflow.data["calculation_uuid"]
    optimized_molecule = rowan.retrieve_calculation_molecules(calculation_uuid)[-1]

    sp_workflow = rowan.submit_workflow(
        initial_molecule=optimized_molecule,
        workflow_type="basic_calculation",
        name=f"{name} {method} single point",
        workflow_data={
            "settings": {
                "method": method,
                "tasks": ["energy"],
                "solvent_settings": {"solvent": "water", "model": "cpcmx"},
            },
        },
    )

    print(f"View workflow privately at: https://labs.rowansci.com/workflow/{sp_workflow.uuid}")
    sp_workflow.wait_for_result().fetch_latest(in_place=True)

    calculation_uuid = sp_workflow.data["calculation_uuid"]
    final_molecule = rowan.retrieve_calculation_molecules(calculation_uuid)[0]
    return final_molecule["energy"]


E1 = compute_energy_with_solvent_correction(
    Molecule.from_smiles("CC(=O)Oc1ccccc1C(=O)O"), Method.OMOL25_CONSERVING_S, "aspirin"
)
print(E1)
