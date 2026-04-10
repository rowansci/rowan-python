"""Descriptors workflow - calculate molecular descriptors."""

from typing import Any

import stjames

from ..folder import Folder
from ..types import MoleculeInput, SolventInput
from ..utils import api_client
from .base import Workflow, WorkflowResult, molecule_to_dict, register_result


@register_result("descriptors")
class DescriptorsResult(WorkflowResult):
    """Result from a descriptors workflow."""

    _stjames_class = stjames.DescriptorsWorkflow

    def __repr__(self) -> str:
        d = self.descriptors or {}
        return f"<DescriptorsResult n={len(d)}>"

    @property
    def descriptors(self) -> dict[str, Any] | None:
        """Computed molecular descriptors."""
        return self._workflow.descriptors


def submit_descriptors_workflow(
    initial_molecule: MoleculeInput,
    solvent: SolventInput = "water",
    do_optimization: bool = True,
    name: str = "Descriptors Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
    is_draft: bool = False,
) -> Workflow:
    """
    Submits a descriptors workflow to the API.

    :param initial_molecule: Molecule to calculate the descriptors of.
    :param solvent: Solvent for COSMO descriptor calculation (e.g. "water").
        When provided, additional COSMO descriptors are computed.
    :param do_optimization: Whether to run GFN2-xTB geometry optimization
        before computing descriptors.
    :param name: Name of the workflow.
    :param folder_uuid: UUID of the folder to place the workflow in.
    :param folder: Folder object to store the workflow in.
    :param max_credits: Maximum number of credits to use for the workflow.
    :param webhook_url: URL that Rowan will POST to when the workflow completes.
    :param is_draft: If True, submit the workflow as a draft without starting execution.
    :returns: Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder:
        folder_uuid = folder.uuid
    mol_dict = molecule_to_dict(initial_molecule)

    workflow_data = {
        "initial_molecule": mol_dict,
        "solvent": solvent,
        "do_optimization": do_optimization,
    }

    workflow = stjames.DescriptorsWorkflow.model_validate(workflow_data)

    data = {
        "workflow_type": "descriptors",
        "workflow_data": workflow.model_dump(serialize_as_any=True, mode="json"),
        "initial_molecule": mol_dict,
        "name": name,
        "folder_uuid": folder_uuid,
        "max_credits": max_credits,
        "webhook_url": webhook_url,
        "is_draft": is_draft,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())
