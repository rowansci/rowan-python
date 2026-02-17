"""MSA workflow - Multiple Sequence Alignment for proteins."""

import stjames

from ..utils import api_client
from .base import Workflow, WorkflowResult, register_result


@register_result("msa")
class MSAResult(WorkflowResult):
    """Result from a Multiple Sequence Alignment (MSA) workflow."""

    _stjames_class = stjames.MSAWorkflow

    def __repr__(self) -> str:
        seqs = getattr(self._workflow, "initial_protein_sequences", []) or []
        return f"<MSAResult sequences={len(seqs)}>"


def submit_msa_workflow(
    initial_protein_sequences: list[stjames.ProteinSequence | str],
    output_formats: list[stjames.MSAFormat],
    name: str = "MSA Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a Multiple Sequence Alignment (MSA) workflow to the API.

    :param initial_protein_sequences: List of protein sequences to align, as ProteinSequence objects
        or strings.
    :param output_formats: List of output formats for the resulting MSA files.
    :param name: The name to assign to the workflow.
    :param folder_uuid: UUID of the folder where the workflow will be stored.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted MSA workflow.
    :raises HTTPError: If the API request fails.
    """
    workflow = stjames.MSAWorkflow(
        initial_protein_sequences=initial_protein_sequences,
        output_formats=output_formats,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "msa",
        "workflow_data": workflow.model_dump(serialize_as_any=True, mode="json"),
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


__all__ = ["MSAResult", "submit_msa_workflow"]
