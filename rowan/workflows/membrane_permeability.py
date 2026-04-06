"""Membrane permeability workflow - predict membrane permeability."""

from typing import Any, Literal

import stjames

from ..folder import Folder
from ..utils import api_client
from .base import (
    MoleculeInput,
    Workflow,
    WorkflowResult,
    molecule_to_dict,
    register_result,
)


@register_result("membrane_permeability")
class MembranePermeabilityResult(WorkflowResult):
    """Result from a membrane-permeability workflow."""

    _stjames_class = stjames.MembranePermeabilityWorkflow

    def __repr__(self) -> str:
        caco = self.caco2_log_p
        bbb = self.bbb_log_p
        return f"<MembranePermeabilityResult caco2_logP={caco} bbb_logP={bbb}>"

    @property
    def caco2_p_app(self) -> float | None:
        """Caco-2 apparent permeability (cm/s)."""
        return getattr(self._workflow, "caco_2_P_app", None)

    @property
    def caco2_log_p(self) -> float | None:
        """Caco-2 log permeability."""
        return getattr(self._workflow, "caco_2_logP", None)

    @property
    def blm_log_p(self) -> float | None:
        """Black lipid membrane log permeability."""
        return getattr(self._workflow, "blm_logP", None)

    @property
    def pampa_log_p(self) -> float | None:
        """PAMPA log permeability."""
        return getattr(self._workflow, "pampa_logP", None)

    @property
    def plasma_log_p(self) -> float | None:
        """Plasma membrane log permeability."""
        return getattr(self._workflow, "plasma_logP", None)

    @property
    def bbb_log_p(self) -> float | None:
        """Blood-brain barrier log permeability."""
        return getattr(self._workflow, "bbb_logP", None)

    @property
    def energy_profile(self) -> list[tuple[float, float]]:
        """Energy profile across the membrane as (position (A), energy (kcal/mol)) pairs."""
        return list(self._workflow.energy_profile)


def submit_membrane_permeability_workflow(
    initial_molecule: MoleculeInput | str,
    method: Literal["gnn-mtl", "pypermm"] = "gnn-mtl",
    name: str = "Membrane Permeability Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
    is_draft: bool = False,
) -> Workflow:
    """
    Submits a membrane-permeability workflow to the API.

    :param initial_molecule: Molecule used in the workflow.
    :param method: Method used to compute membrane permeability.
    :param name: Name of the workflow.
    :param folder_uuid: UUID of the folder to store the workflow in.
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
    data: dict[str, Any] = {
        "workflow_type": "membrane_permeability",
        "name": name,
        "folder_uuid": folder_uuid,
        "max_credits": max_credits,
        "webhook_url": webhook_url,
        "is_draft": is_draft,
    }

    match method:
        case "gnn-mtl":
            if not isinstance(initial_molecule, str):
                raise ValueError("initial_molecule must be a SMILES string for the gnn-mtl method")
            workflow = stjames.MembranePermeabilityWorkflow(
                initial_smiles=initial_molecule,
                membrane_permeability_method="chemprop_ohlsson2025",
            )

            data["initial_smiles"] = initial_molecule

        case "pypermm":
            if isinstance(initial_molecule, str):
                raise ValueError("Cannot specify molecule as SMILES for PyPermm")
            initial_molecule = molecule_to_dict(initial_molecule)

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
