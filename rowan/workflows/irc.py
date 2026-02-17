"""IRC workflow - Intrinsic Reaction Coordinate calculations."""

from typing import Any

import stjames
from rdkit import Chem

from ..utils import api_client
from .base import RdkitMol, StJamesMolecule, Workflow, WorkflowResult, register_result


@register_result("irc")
class IRCResult(WorkflowResult):
    """Result from an Intrinsic Reaction Coordinate (IRC) workflow."""

    _stjames_class = stjames.IRCWorkflow

    def __repr__(self) -> str:
        n_fwd = len(self.irc_forward)
        n_bwd = len(self.irc_backward)
        return f"<IRCResult forward_steps={n_fwd} backward_steps={n_bwd}>"

    @property
    def starting_ts(self) -> str | None:
        """UUID of optimized TS before IRC."""
        return self._workflow.starting_TS

    @property
    def irc_forward(self) -> list[str]:
        """Forward IRC path UUIDs."""
        return list(self._workflow.irc_forward)

    @property
    def irc_backward(self) -> list[str]:
        """Backward IRC path UUIDs."""
        return list(self._workflow.irc_backward)


def submit_irc_workflow(
    initial_molecule: dict[str, Any] | StJamesMolecule | RdkitMol | None = None,
    method: stjames.Method | str = "uma_m_omol",
    preopt: bool = True,
    step_size: float = 0.05,
    max_irc_steps: int = 30,
    name: str = "IRC Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits an Intrinsic Reaction Coordinate (IRC) workflow to the API.

    :param initial_molecule: The initial molecule to perform the IRC calculation on.
    :param method: The computational method to use for the IRC calculation.
    :param preopt: Whether to perform a pre-optimization of the molecule.
    :param step_size: The step size to use for the IRC calculation.
    :param max_irc_steps: The maximum number of IRC steps to perform.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted IRC workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if isinstance(initial_molecule, StJamesMolecule):
        initial_molecule = initial_molecule.model_dump(mode="json")
    elif isinstance(initial_molecule, Chem.rdchem.Mol | Chem.rdchem.RWMol):
        initial_molecule = StJamesMolecule.from_rdkit(initial_molecule, cid=0)

    if isinstance(method, str):
        method = stjames.Method(method)

    workflow = stjames.IRCWorkflow(
        initial_molecule=initial_molecule,
        settings=stjames.Settings(
            method=method,
            tasks=[],
            corrections=[],
            mode="auto",
        ),
        preopt=preopt,
        step_size=step_size,
        max_irc_steps=max_irc_steps,
        mode="manual",
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "irc",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": initial_molecule,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


__all__ = ["IRCResult", "submit_irc_workflow"]
