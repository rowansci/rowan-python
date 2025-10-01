# ruff: noqa
from stjames import Molecule, Method

import rowan

# rowan.api_key = ""


def compute_energy_with_solvent_correction(molecule: Molecule, method: Method, name: str) -> float:
    opt = rowan.submit_workflow(
        initial_molecule=molecule,
        workflow_type="basic_calculation",
        name=f"{name} {method} optimization",
        workflow_data={
            "settings": {"method": method, "tasks": ["optimize"]},
        },
    )

    opt.wait_for_result()
    opt.fetch_latest(in_place=True)

    calculation_uuid = opt.data["calculation_uuid"]
    optimized_molecule = rowan.retrieve_calculation_molecules(calculation_uuid)[-1]

    sp = rowan.submit_workflow(
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

    sp.wait_for_result()
    sp.fetch_latest(in_place=True)

    calculation_uuid = sp.data["calculation_uuid"]
    final_molecule = rowan.retrieve_calculation_molecules(calculation_uuid)[0]
    return final_molecule["energy"]


E1 = compute_energy_with_solvent_correction(
    Molecule.from_smiles("CC(=O)Oc1ccccc1C(=O)O"), Method.OMOL25_CONSERVING_S, "aspirin"
)
print(E1)
