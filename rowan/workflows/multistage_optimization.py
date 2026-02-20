"""Multistage optimization workflow - optimize molecules with staged methods."""

from typing import Any

import stjames
from rdkit import Chem

from ..utils import api_client
from .base import RdkitMol, StJamesMolecule, Workflow, WorkflowResult, register_result


@register_result("multistage_opt")
class MultiStageOptResult(WorkflowResult):
    """Result from a multistage optimization workflow."""

    _stjames_class = stjames.MultiStageOptWorkflow

    def __repr__(self) -> str:
        final_energy = getattr(self._workflow, "final_energy", None)
        return f"<MultiStageOptResult final_energy={final_energy}>"


def submit_multistage_optimization_workflow(
    initial_molecule: dict[str, Any] | stjames.Molecule | RdkitMol,
    mode: str = "rapid",
    solvent: str | None = None,
    xtb_preopt: bool = True,
    transition_state: bool = False,
    frequencies: bool = False,
    name: str = "Multistage Optimization Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a multistage optimization workflow to the API.

    :param initial_molecule: The molecule to optimize.
    :param mode: The mode to run the calculation in.
    :param solvent: The solvent to use for the final single-point calculation.
    :param xtb_preopt: Whether to pre-optimize with xTB.
    :param transition_state: Whether this is a transition state optimization.
    :param frequencies: Whether to calculate frequencies.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if isinstance(initial_molecule, StJamesMolecule):
        initial_molecule = initial_molecule.model_dump(mode="json")
    elif isinstance(initial_molecule, Chem.rdchem.Mol | Chem.rdchem.RWMol):
        initial_molecule = stjames.Molecule.from_rdkit(initial_molecule, cid=0).model_dump(
            mode="json"
        )

    workflow = stjames.MultiStageOptWorkflow(
        initial_molecule=initial_molecule,
        mode=mode,
        solvent=solvent,
        xtb_preopt=xtb_preopt,
        transition_state=transition_state,
        frequencies=frequencies,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "multistage_opt",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": initial_molecule,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


__all__ = ["MultiStageOptResult", "submit_multistage_optimization_workflow"]
