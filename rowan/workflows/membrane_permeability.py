"""Membrane permeability workflow - predict membrane permeability."""

from typing import Any, Literal

import stjames
from rdkit import Chem

from ..utils import api_client
from .base import RdkitMol, StJamesMolecule, Workflow, WorkflowResult, register_result


@register_result("membrane_permeability")
class MembranePermeabilityResult(WorkflowResult):
    """Result from a membrane permeability workflow."""

    _stjames_class = stjames.MembranePermeabilityWorkflow

    def __repr__(self) -> str:
        perm = getattr(self._workflow, "permeability", None)
        return f"<MembranePermeabilityResult permeability={perm}>"


def submit_membrane_permeability_workflow(
    initial_molecule: dict[str, Any] | StJamesMolecule | RdkitMol | str,
    method: Literal["gnn-mtl", "pypermm"] = "gnn-mtl",
    name: str = "Membrane Permeability Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a membrane permeability workflow to the API.

    :param initial_molecule: The molecule used in the workflow.
    :param method: The method used to compute membrane permeability.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to store the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    data: dict[str, Any] = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "membrane_permeability",
        "max_credits": max_credits,
    }

    match method:
        case "gnn-mtl":
            assert isinstance(initial_molecule, str)
            workflow = stjames.MembranePermeabilityWorkflow(
                initial_smiles=initial_molecule,
                membrane_permeability_method="chemprop_ohlsson2025",
            )

            data["initial_smiles"] = initial_molecule

        case "pypermm":
            if isinstance(initial_molecule, str):
                raise ValueError("Cannot specify molecule as SMILES for PyPermm")
            elif isinstance(initial_molecule, StJamesMolecule):
                initial_molecule = initial_molecule.model_dump(mode="json")
            elif isinstance(initial_molecule, Chem.rdchem.Mol | Chem.rdchem.RWMol):
                initial_molecule = StJamesMolecule.from_rdkit(initial_molecule, cid=0).model_dump(
                    mode="json"
                )

            workflow = stjames.MembranePermeabilityWorkflow(
                initial_molecule=initial_molecule,
                membrane_permeability_method="pypermm",
            )

            data["initial_molecule"] = initial_molecule

        case _:
            raise ValueError(f"Unexpected {method=}")

    data["workflow_data"] = workflow.model_dump(mode="json")

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


__all__ = ["MembranePermeabilityResult", "submit_membrane_permeability_workflow"]
