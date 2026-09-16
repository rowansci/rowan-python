"""ADMET workflow - Absorption, Distribution, Metabolism, Excretion, and Toxicity."""

import stjames

from ..folder import Folder
from ..utils import api_client
from .base import SMILES, Workflow, WorkflowResult, extract_smiles, register_result


@register_result("admet")
class ADMETResult(WorkflowResult):
    """Result from an ADMET workflow."""

    _stjames_class = stjames.ADMETWorkflow

    def __repr__(self) -> str:
        props = self.properties or {}
        preview = {k: props[k] for k in list(props.keys())[:5]} if props else {}
        return f"<ADMETResult properties={len(props)} preview={preview}>"

    @property
    def properties(self) -> dict[str, float | int] | None:
        """ADMET properties (molecular weight, logP, TPSA, etc.)."""
        return self._workflow.properties


def submit_admet_workflow(
    initial_smiles: SMILES,
    name: str = "ADMET Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
    is_draft: bool = False,
) -> Workflow[ADMETResult]:
    """Submits an ADMET workflow to predict drug-likeness properties.

    Args:
        initial_smiles: molecule with SMILES for ADMET prediction
        name: name of the workflow
        folder_uuid: UUID of the folder to store the workflow in
        folder: destination folder
        max_credits: maximum credits for the workflow
        webhook_url: URL that Rowan will POST to when the workflow completes
        is_draft: save as a draft without starting execution

    Returns:
        submitted workflow

    Raises:
        ValueError: molecule has no SMILES associated with it
        httpx.HTTPStatusError: request to the API fails
    """
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder:
        folder_uuid = folder.uuid
    initial_smiles = extract_smiles(initial_smiles)
    workflow = stjames.ADMETWorkflow(initial_smiles=initial_smiles)

    data = {
        "workflow_type": "admet",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_smiles": initial_smiles,
        "name": name,
        "folder_uuid": folder_uuid,
        "max_credits": max_credits,
        "webhook_url": webhook_url,
        "is_draft": is_draft,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow[ADMETResult](**response.json())
