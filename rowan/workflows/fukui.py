"""Fukui workflow - calculate Fukui indices for reactivity prediction."""

import stjames

from ..folder import Folder
from ..utils import api_client
from .base import (
    StructureInput,
    Workflow,
    WorkflowResult,
    molecule_to_dict,
    register_result,
    require_coordinates,
)


@register_result("fukui")
class FukuiResult(WorkflowResult):
    """Result from a Fukui index workflow."""

    _stjames_class = stjames.FukuiIndexWorkflow

    def __repr__(self) -> str:
        gei = self.global_electrophilicity_index
        return f"<FukuiResult global_electrophilicity_index={gei} eV>"

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
    initial_molecule: StructureInput,
    optimization_method: str = "gfn2_xtb",
    fukui_method: str = "gfn1_xtb",
    solvent_settings: dict[str, str] | None = None,
    name: str = "Fukui Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
    is_draft: bool = False,
) -> Workflow[FukuiResult]:
    """Submits a Fukui workflow to the API.

    Args:
        initial_molecule: molecule to calculate the Fukui indices of
        optimization_method: method to use for the optimization
        fukui_method: method to use for the Fukui calculation
        solvent_settings: optional implicit solvent for the Fukui calculation. A dict
            with two keys:

            - `"solvent"`: solvent name string (e.g. `"water"`, `"dichloromethane"`,
              `"dmso"`). See `rowan.Solvent` for all valid values.
            - `"model"`: solvation model (e.g. `"alpb"`, `"gbsa"`, `"cpcmx"` for xTB;
              `"cpcm"`, `"pcm"` for DFT). Must be compatible with the engine for the chosen method

            Example: `solvent_settings={"solvent": "water", "model": "alpb"}`
        name: name of the workflow
        folder_uuid: UUID of the folder to place the workflow in
        folder: destination folder
        max_credits: maximum credits for the workflow
        webhook_url: URL that Rowan will POST to when the workflow completes
        is_draft: save as a draft without starting execution

    Returns:
        submitted workflow

    Raises:
        httpx.HTTPStatusError: request to the API fails
    """
    require_coordinates(initial_molecule)
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder:
        folder_uuid = folder.uuid
    mol_dict = molecule_to_dict(initial_molecule)

    optimization_settings = stjames.Settings(method=optimization_method)
    fukui_settings = stjames.Settings(method=fukui_method, solvent_settings=solvent_settings)

    workflow_data = {
        "opt_settings": optimization_settings.model_dump(mode="json"),
        "opt_engine": stjames.Method(optimization_method).default_engine(),
        "fukui_settings": fukui_settings.model_dump(mode="json"),
        "fukui_engine": stjames.Method(fukui_method).default_engine(),
    }

    data = {
        "workflow_type": "fukui",
        "workflow_data": workflow_data,
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
        return Workflow[FukuiResult](**response.json())
