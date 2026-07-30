"""logP workflow - predict the octanol/water partition coefficient."""

from typing import Literal

import stjames

from ..folder import Folder
from ..utils import api_client
from .base import SMILES, Workflow, WorkflowResult, extract_smiles, register_result


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
    initial_smiles: SMILES,
    method: Literal["chemprop_sangster2026", "crippen", "cosmors"] = "chemprop_sangster2026",
    name: str = "LogP Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
    is_draft: bool = False,
) -> Workflow:
    """
    Submits a logP workflow to the API.

    :param initial_smiles: Molecule to predict logP for. Accepts a SMILES string or any
        molecule type (RowanMolecule, stjames.Molecule, RDKit Mol, or dict). The molecule
        must have a SMILES string associated with it, as this workflow is SMILES-based and
        does not use 3D coordinates.
    :param method: logP prediction method:
        - "chemprop_sangster2026": chemprop v2 D-MPNN trained on experimental octanol/water
          logP from the Sangster dataset.
        - "crippen": RDKit implementation of the Wildman-Crippen atom-contribution model.
        - "cosmors": Boltzmann-weighted COSMO-RS logP over a conformer ensemble.
    :param name: Name of the workflow.
    :param folder_uuid: UUID of the folder to store the workflow in.
    :param folder: Folder object to store the workflow in.
    :param max_credits: Maximum number of credits to use for the workflow.
    :param webhook_url: URL that Rowan will POST to when the workflow completes.
    :param is_draft: If True, submit the workflow as a draft without starting execution.
    :returns: Workflow object representing the submitted workflow.
    :raises ValueError: If the molecule has no SMILES associated with it.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder:
        folder_uuid = folder.uuid
    initial_smiles = extract_smiles(initial_smiles)
    workflow = stjames.LogPWorkflow(initial_smiles=initial_smiles, logp_method=method)

    data = {
        "workflow_type": "logp",
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
        return Workflow(**response.json())
