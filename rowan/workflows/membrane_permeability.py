"""Membrane permeability workflow - predict membrane permeability."""

from typing import Any, Literal

import stjames

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
    """Result from a membrane permeability workflow."""

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
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a membrane permeability workflow to the API.

    :param initial_molecule: Molecule used in the workflow.
    :param method: Method used to compute membrane permeability.
    :param name: Name of the workflow.
    :param folder_uuid: UUID of the folder to store the workflow in.
    :param max_credits: Maximum number of credits to use for the workflow.
    :returns: Workflow object representing the submitted workflow.
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
