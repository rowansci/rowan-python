"""logP workflow - predict the octanol/water partition coefficient."""

from typing import Any, Literal

import stjames
from rowan.folder import Folder
from rowan.types import SMILES, StructureInput
from rowan.utils import api_client

from .base import (
    Workflow,
    WorkflowResult,
    extract_smiles,
    molecule_to_stjames,
    register_result,
    require_coordinates,
)


@register_result("logp")
class LogPResult(WorkflowResult):
    """Result from a logP workflow."""

    _stjames_class = stjames.LogPWorkflow

    def __repr__(self) -> str:
        return f"<LogPResult logp={self.logp}>"

    @property
    def logp(self) -> float | None:
        """Predicted octanol/water partition coefficient."""
        return self._workflow.logp


def submit_logp_workflow(
    initial_smiles: SMILES | StructureInput | dict[str, Any] | None = None,
    method: Literal["chemprop_sangster2026", "crippen", "cosmors"] = "chemprop_sangster2026",
    name: str = "LogP Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
    is_draft: bool = False,
    *,
    initial_molecule: StructureInput | dict[str, Any] | None = None,
) -> Workflow[LogPResult]:
    """Submits a logP workflow to the API.

    Args:
        initial_smiles: SMILES string or molecule with associated SMILES
        initial_molecule: nonperiodic 3D molecule to optimize and score as a single conformer
            with COSMO-RS. Provide exactly one of `initial_smiles` and `initial_molecule`
        method: logP prediction method:
            - "chemprop_sangster2026": chemprop v2 D-MPNN trained on experimental octanol/water
              logP from the Sangster dataset.
            - "crippen": RDKit implementation of the Wildman-Crippen atom-contribution model.
            - "cosmors": Boltzmann-weighted COSMO-RS logP over a conformer ensemble
        name: name of the workflow
        folder_uuid: UUID of the folder to store the workflow in
        folder: destination folder
        max_credits: maximum credits for the workflow
        webhook_url: URL that Rowan will POST to when the workflow completes
        is_draft: save as a draft without starting execution

    Returns:
        submitted workflow

    Raises:
        ValueError: molecular input is missing, ambiguous, or incompatible with the method
        httpx.HTTPStatusError: request to the API fails
    """
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder:
        folder_uuid = folder.uuid
    smiles = extract_smiles(initial_smiles) if initial_smiles is not None else ""
    molecule = None
    if initial_molecule is not None:
        if not isinstance(initial_molecule, dict):
            require_coordinates(initial_molecule)
        molecule = molecule_to_stjames(initial_molecule)
    workflow = stjames.LogPWorkflow(
        initial_smiles=smiles, initial_molecule=molecule, logp_method=method
    )
    workflow_data = workflow.model_dump(mode="json")
    data = {
        "workflow_type": "logp",
        "workflow_data": workflow_data,
        "name": name,
        "folder_uuid": folder_uuid,
        "max_credits": max_credits,
        "webhook_url": webhook_url,
        "is_draft": is_draft,
    }
    if molecule is not None:
        data["initial_molecule"] = workflow_data["initial_molecule"]
    else:
        data["initial_smiles"] = smiles

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow[LogPResult](**response.json())
