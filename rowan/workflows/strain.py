"""Strain workflow - calculate molecular strain energy."""

from typing import Any

import stjames
from rdkit import Chem

from ..utils import api_client
from .base import RdkitMol, StJamesMolecule, Workflow, WorkflowResult, register_result


@register_result("strain")
class StrainResult(WorkflowResult):
    """Result from a strain workflow."""

    _stjames_class = stjames.StrainWorkflow

    def __repr__(self) -> str:
        return f"<StrainResult strain={self.strain}>"

    @property
    def strain(self) -> float | None:
        """Computed strain energy (kcal/mol)."""
        return self._workflow.strain

    @property
    def conformers(self) -> list[str | None]:
        """UUIDs of conformers."""
        return list(self._workflow.conformers)


def submit_strain_workflow(
    initial_molecule: dict[str, Any] | stjames.Molecule | RdkitMol,
    name: str = "Strain Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a strain workflow to the API.

    :param initial_molecule: The molecule used in the scan.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to store the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if isinstance(initial_molecule, StJamesMolecule):
        initial_molecule = initial_molecule.model_dump(mode="json")
    elif isinstance(initial_molecule, Chem.rdchem.Mol | Chem.rdchem.RWMol):
        initial_molecule = stjames.Molecule.from_rdkit(initial_molecule, cid=0)

    workflow = stjames.StrainWorkflow(initial_molecule=initial_molecule)

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "strain",
        "workflow_data": workflow.model_dump(serialize_as_any=True, mode="json"),
        "initial_molecule": initial_molecule,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


__all__ = ["StrainResult", "submit_strain_workflow"]
