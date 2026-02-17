"""Descriptors workflow - calculate molecular descriptors."""

from typing import Any

import stjames
from rdkit import Chem

from ..utils import api_client
from .base import RdkitMol, StJamesMolecule, Workflow, WorkflowResult, register_result


@register_result("descriptors")
class DescriptorsResult(WorkflowResult):
    """Result from a descriptors workflow."""

    _stjames_class = stjames.DescriptorsWorkflow

    def __repr__(self) -> str:
        d = self.descriptors or {}
        preview = {k: d[k] for k in list(d.keys())[:5]} if d else {}
        return f"<DescriptorsResult count={len(d)} preview={preview}>"

    @property
    def descriptors(self) -> dict | None:
        """Computed molecular descriptors."""
        return self._workflow.descriptors


def submit_descriptors_workflow(
    initial_molecule: dict[str, Any] | StJamesMolecule | RdkitMol,
    name: str = "Descriptors Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a descriptors workflow to the API.

    :param initial_molecule: The molecule to calculate the descriptors of.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if isinstance(initial_molecule, StJamesMolecule):
        initial_molecule = initial_molecule.model_dump(mode="json")
    elif isinstance(initial_molecule, Chem.rdchem.Mol | Chem.rdchem.RWMol):
        initial_molecule = StJamesMolecule.from_rdkit(initial_molecule, cid=0)

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


__all__ = ["DescriptorsResult", "submit_descriptors_workflow"]
