"""Redox potential workflow - calculate oxidation/reduction potentials."""

from typing import Any

import stjames
from rdkit import Chem

from ..utils import api_client
from .base import RdkitMol, StJamesMolecule, Workflow, WorkflowResult, register_result


@register_result("redox_potential")
class RedoxPotentialResult(WorkflowResult):
    """Result from a redox potential workflow."""

    _stjames_class = stjames.RedoxPotentialWorkflow

    def __repr__(self) -> str:
        ox = self.oxidation_potential
        red = self.reduction_potential
        return f"<RedoxPotentialResult oxidation={ox} reduction={red}>"

    @property
    def oxidation_potential(self) -> float | None:
        """Oxidation potential in V."""
        return self._workflow.oxidation_potential

    @property
    def reduction_potential(self) -> float | None:
        """Reduction potential in V."""
        return self._workflow.reduction_potential


def submit_redox_potential_workflow(
    initial_molecule: dict[str, Any] | stjames.Molecule | RdkitMol,
    reduction: bool = False,
    oxidization: bool = True,
    mode: str = "rapid",
    name: str = "Redox Potential Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a redox potential workflow to the API.

    :param initial_molecule: The molecule to calculate the redox potential of.
    :param reduction: Whether to calculate the reduction potential.
    :param oxidization: Whether to calculate the oxidization potential.
    :param mode: The mode to run the calculation in.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if isinstance(initial_molecule, StJamesMolecule):
        initial_molecule = initial_molecule.model_dump(mode="json")
    elif isinstance(initial_molecule, Chem.rdchem.Mol | Chem.rdchem.RWMol):
        initial_molecule = stjames.Molecule.from_rdkit(initial_molecule, cid=0)

    workflow = stjames.RedoxPotentialWorkflow(
        initial_molecule=initial_molecule,
        oxidation=oxidization,
        reduction=reduction,
        mode=mode,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "redox_potential",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": initial_molecule,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


__all__ = ["RedoxPotentialResult", "submit_redox_potential_workflow"]
