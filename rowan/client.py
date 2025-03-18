import time
from typing import Any

import stjames

from .workflow import Workflow

""" A high-level interface to submitting a calculation. """


def compute(
    molecule: str | stjames.Molecule,
    workflow_type: str,
    name: str = "",
    folder_uuid: stjames.UUID | None = None,
    blocking: bool = True,
    ping_interval: int = 5,
    **workflow_data,
) -> dict[str, Any]:
    """
    High-level function to compute and return workflows.

    :param molecule: Molecule to compute
    :param workflow_type: name of workflow to compute
    :param name: name for the job
    :param folder_uuid: folder to store the job
    :param blocking: whether to wait for the job to finish
    :param ping_interval: interval to check if the job is finished
    :param workflow_data: additional data to pass to the workflow
    :return: result of the workflow
    """

    if isinstance(molecule, str):
        stjmol = stjames.Molecule.from_smiles(molecule)
    elif isinstance(molecule, stjames.Molecule):
        stjmol = molecule
    else:
        raise ValueError("Invalid type for `molecule`!")

    result = Workflow.submit(
        initial_molecule=stjmol,
        workflow_type=workflow_type,
        name=name,
        folder_uuid=folder_uuid,
        workflow_data=workflow_data,
    )

    if blocking:
        uuid = result["uuid"]
        while not Workflow.is_finished(uuid):
            time.sleep(ping_interval)

        return Workflow.retrieve(uuid)

    return result
