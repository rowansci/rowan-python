"""ADMET workflow - Absorption, Distribution, Metabolism, Excretion, and Toxicity."""

import stjames

from ..folder import Folder
from ..types import MoleculeInput
from ..utils import api_client
from .base import Workflow, WorkflowResult, extract_smiles, register_result


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
    initial_smiles: str | MoleculeInput,
    name: str = "ADMET Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
) -> Workflow:
    """
    Submits an ADMET workflow to predict drug-likeness properties.

    :param initial_smiles: Molecule to calculate ADMET properties for. Accepts a
        SMILES string or any molecule type (RowanMolecule, stjames.Molecule, RDKit Mol,
        or dict). The molecule must have a SMILES string associated with it, as ADMET
        models are 2D/SMILES-based and do not use 3D coordinates.
    :param name: Name of the workflow.
    :param folder_uuid: UUID of the folder to store the workflow in.
    :param folder: Folder object to store the workflow in.
    :param max_credits: Maximum number of credits to use for the workflow.
    :param webhook_url: URL that Rowan will POST to when the workflow completes.
    :returns: Workflow object representing the submitted workflow.
    :raises ValueError: If the molecule has no SMILES associated with it.
    :raises requests.HTTPError: if the request to the API fails.
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
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())
