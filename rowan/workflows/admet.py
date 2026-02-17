"""ADMET workflow - Absorption, Distribution, Metabolism, Excretion, and Toxicity."""

import stjames

from ..utils import api_client
from .base import Workflow, WorkflowResult, register_result


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
    initial_smiles: str,
    name: str = "ADMET Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits an ADMET workflow to predict drug-likeness properties.

    :param initial_smiles: The molecule SMILES to calculate ADMET properties for.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to store the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    workflow = stjames.ADMETWorkflow(initial_smiles=initial_smiles)

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "admet",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_smiles": initial_smiles,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


__all__ = ["ADMETResult", "submit_admet_workflow"]
