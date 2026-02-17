"""Hydrogen bond basicity workflow - calculate hydrogen bond basicity."""

from typing import Any

import stjames
from rdkit import Chem

from ..utils import api_client
from .base import RdkitMol, StJamesMolecule, Workflow, WorkflowResult, register_result


@register_result("hydrogen_bond_basicity")
class HydrogenBondBasicityResult(WorkflowResult):
    """Result from a hydrogen bond basicity workflow."""

    _stjames_class = stjames.HydrogenBondBasicityWorkflow

    def __repr__(self) -> str:
        basicity = getattr(self._workflow, "hydrogen_bond_basicity", None)
        return f"<HydrogenBondBasicityResult basicity={basicity}>"


def submit_hydrogen_bond_basicity_workflow(
    initial_molecule: dict[str, Any] | StJamesMolecule | RdkitMol,
    do_csearch: bool = True,
    do_optimization: bool = True,
    name: str = "Hydrogen Bond Basicity Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a hydrogen bond basicity workflow to the API.

    :param initial_molecule: The molecule to calculate hydrogen bond basicity for.
    :param do_csearch: Whether to perform a conformational search.
    :param do_optimization: Whether to perform an optimization.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if isinstance(initial_molecule, StJamesMolecule):
        initial_molecule = initial_molecule.model_dump(mode="json")
    elif isinstance(initial_molecule, Chem.rdchem.Mol | Chem.rdchem.RWMol):
        initial_molecule = StJamesMolecule.from_rdkit(initial_molecule, cid=0).model_dump(
            mode="json"
        )

    workflow = stjames.HydrogenBondBasicityWorkflow(
        initial_molecule=initial_molecule,
        do_csearch=do_csearch,
        do_optimization=do_optimization,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "hydrogen_bond_basicity",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": initial_molecule,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


__all__ = ["HydrogenBondBasicityResult", "submit_hydrogen_bond_basicity_workflow"]
