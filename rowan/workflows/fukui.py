"""Fukui workflow - calculate Fukui indices for reactivity prediction."""

import stjames

from ..folder import Folder
from ..types import MoleculeInput
from ..utils import api_client
from .base import Workflow, WorkflowResult, molecule_to_dict, register_result


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
    initial_molecule: MoleculeInput,
    optimization_method: str = "gfn2_xtb",
    fukui_method: str = "gfn1_xtb",
    solvent_settings: dict[str, str] | None = None,
    name: str = "Fukui Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
) -> Workflow:
    """
    Submits a Fukui workflow to the API.

    :param initial_molecule: Molecule to calculate the Fukui indices of.
    :param optimization_method: Method to use for the optimization.
    :param fukui_method: Method to use for the Fukui calculation.
    :param solvent_settings: Optional implicit solvent for the Fukui calculation. A dict
        with two keys:

        - ``"solvent"``: solvent name string (e.g. ``"water"``, ``"dichloromethane"``,
          ``"dmso"``). See ``rowan.Solvent`` for all valid values.
        - ``"model"``: solvation model. Use ``"alpb"`` or ``"gbsa"`` for xTB methods
          (the defaults ``gfn1_xtb`` / ``gfn2_xtb``); use ``"cpcm"`` or ``"smd"`` for
          DFT methods.

        Example: ``solvent_settings={"solvent": "water", "model": "alpb"}``
    :param name: Name of the workflow.
    :param folder_uuid: UUID of the folder to place the workflow in.
    :param folder: Folder object to store the workflow in.
    :param max_credits: Maximum number of credits to use for the workflow.
    :param webhook_url: URL that Rowan will POST to when the workflow completes.
    :returns: Workflow object representing the submitted workflow.
    :raises ValueError: If the solvent model is incompatible with the chosen method.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder:
        folder_uuid = folder.uuid
    if solvent_settings is not None:
        model = solvent_settings.get("model")
        is_xtb = stjames.Method(fukui_method) in stjames.XTB_METHODS
        xtb_models = {"alpb", "gbsa"}
        dft_models = {"pcm", "cpcm", "cosmo", "cpcmx", "smd"}
        if is_xtb and model not in xtb_models:
            raise ValueError(
                f"xTB Fukui methods require 'alpb' or 'gbsa' solvation model, got '{model}'"
            )
        if not is_xtb and model not in dft_models:
            raise ValueError(
                f"DFT Fukui methods require 'cpcm', 'smd', 'pcm', 'cosmo', or 'cpcmx' "
                f"solvation model, got '{model}'"
            )

    initial_molecule = molecule_to_dict(initial_molecule)

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
        "initial_molecule": initial_molecule,
        "name": name,
        "folder_uuid": folder_uuid,
        "max_credits": max_credits,
        "webhook_url": webhook_url,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())
