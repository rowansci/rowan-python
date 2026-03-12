"""Descriptors workflow - calculate molecular descriptors."""

from typing import Any

import stjames

from ..types import MoleculeInput
from ..utils import api_client
from .base import Workflow, WorkflowResult, molecule_to_dict, register_result


@register_result("descriptors")
class DescriptorsResult(WorkflowResult):
    """Result from a descriptors workflow."""

    _stjames_class = stjames.DescriptorsWorkflow

    def __repr__(self) -> str:
        d = self.descriptors or {}
        preview = {k: d[k] for k in list(d.keys())[:5]} if d else {}
        return f"<DescriptorsResult count={len(d)} preview={preview}>"

    @property
    def descriptors(self) -> dict[str, Any] | None:
        """Computed molecular descriptors."""
        return self._workflow.descriptors


def submit_descriptors_workflow(
    initial_molecule: MoleculeInput,
    name: str = "Descriptors Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a descriptors workflow to the API.

    :param initial_molecule: Molecule to calculate the descriptors of.
    :param name: Name of the workflow.
    :param folder_uuid: UUID of the folder to place the workflow in.
    :param max_credits: Maximum number of credits to use for the workflow.
    :returns: Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    initial_molecule = molecule_to_dict(initial_molecule)

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "descriptors",
        "workflow_data": {},
        "initial_molecule": initial_molecule,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())
