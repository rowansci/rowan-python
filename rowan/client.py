import cctk
import stjames
import time
from typing import Optional

from .utils import cctk_to_stjames, smiles_to_stjames
from .workflow import Workflow

""" A high-level interface to submitting a calculation. """


def compute(
    molecule: str | cctk.Molecule | stjames.Molecule,
    workflow_type: str,
    name: str = "",
    folder_uuid: Optional[stjames.UUID] = None,
    blocking: bool = True,
    ping_interval: int = 5,
    **workflow_data,
) -> dict:
    """High-level function to compute and return workflows."""

    if isinstance(molecule, str):
        stjmol = smiles_to_stjames(molecule)
    elif isinstance(molecule, cctk.Molecule):
        stjmol = cctk_to_stjames(molecule)
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

        completed_result = Workflow.retrieve(uuid)
        return completed_result

    else:
        return result
