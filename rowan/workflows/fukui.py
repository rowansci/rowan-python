"""Fukui workflow - calculate Fukui indices for reactivity prediction."""

from typing import Any

import stjames
from rdkit import Chem

from ..utils import api_client
from .base import RdkitMol, StJamesMolecule, Workflow, WorkflowResult, register_result


@register_result("fukui")
class FukuiResult(WorkflowResult):
    """Result from a Fukui index workflow."""

    _stjames_class = stjames.FukuiIndexWorkflow

    def __repr__(self) -> str:
        gei = self.global_electrophilicity_index
        return f"<FukuiResult global_electrophilicity_index={gei}>"

    @property
    def global_electrophilicity_index(self) -> float | None:
        """Global electrophilicity index."""
        return self._workflow.global_electrophilicity_index

    @property
    def fukui_positive(self) -> list[float] | None:
        """Fukui f+ indices (electrophilic attack susceptibility)."""
        return list(self._workflow.fukui_positive) if self._workflow.fukui_positive else None

    @property
    def fukui_negative(self) -> list[float] | None:
        """Fukui f- indices (nucleophilic attack susceptibility)."""
        return list(self._workflow.fukui_negative) if self._workflow.fukui_negative else None

    @property
    def fukui_zero(self) -> list[float] | None:
        """Fukui f0 indices (radical attack susceptibility)."""
        return list(self._workflow.fukui_zero) if self._workflow.fukui_zero else None


def submit_fukui_workflow(
    initial_molecule: dict[str, Any] | StJamesMolecule | RdkitMol,
    optimization_method: str = "gfn2_xtb",
    fukui_method: str = "gfn1_xtb",
    solvent_settings: dict[str, str] | None = None,
    name: str = "Fukui Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a fukui workflow to the API.

    :param initial_molecule: The molecule to calculate the fukui indices of.
    :param optimization_method: The method to use for the optimization.
    :param fukui_method: The method to use for the fukui calculation.
    :param solvent_settings: The solvent settings to use for the fukui calculation.
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

    optimization_settings = stjames.Settings(method=optimization_method)
    fukui_settings = stjames.Settings(method=fukui_method, solvent_settings=solvent_settings)

    stjames.FukuiIndexWorkflow(
        initial_molecule=initial_molecule,
        optimization_settings=optimization_settings,
        fukui_settings=fukui_settings,
    )

    workflow_data = {
        "opt_settings": optimization_settings.model_dump(mode="json"),
        "opt_engine": stjames.Method(optimization_method).default_engine(),
        "fukui_settings": fukui_settings.model_dump(mode="json"),
        "fukui_engine": stjames.Method(fukui_method).default_engine(),
    }

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "fukui",
        "workflow_data": workflow_data,
        "initial_molecule": initial_molecule,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


__all__ = ["FukuiResult", "submit_fukui_workflow"]
