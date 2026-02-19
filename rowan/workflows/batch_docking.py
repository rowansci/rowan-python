"""Batch docking workflow - high-throughput molecular docking."""

import stjames

from ..protein import Protein
from ..utils import api_client
from .base import Workflow, WorkflowResult, register_result


@register_result("batch_docking")
class BatchDockingResult(WorkflowResult):
    """Result from a batch docking workflow."""

    _stjames_class = stjames.BatchDockingWorkflow

    def __repr__(self) -> str:
        scores = self.scores
        valid = {k: v for k, v in scores.items() if v is not None}
        if valid:
            best_smiles = min(valid.keys(), key=lambda k: valid[k])
            best_score = valid[best_smiles]
        else:
            best_smiles, best_score = None, None
        return f"<BatchDockingResult ligands={len(scores)} best=({best_score}, {best_smiles!r})>"

    @property
    def scores(self) -> dict[str, float | None]:
        """Docking scores indexed by SMILES."""
        smiles_list = self._workflow.initial_smiles_list
        scores = self._workflow.best_scores
        return dict(zip(smiles_list, scores, strict=True))


def submit_batch_docking_workflow(
    smiles_list: list[str],
    protein: str | Protein,
    pocket: list[list[float]],
    executable: str = "qvina2",
    scoring_function: str = "vina",
    exhaustiveness: float = 8,
    name: str = "Batch Docking Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a batch docking workflow to the API.

    :param smiles_list: The SMILES strings to dock.
    :param protein: The protein to dock (UUID or Protein object).
    :param pocket: Binding pocket coordinates [[x,y,z], [x,y,z]].
    :param executable: Which docking implementation to use.
    :param scoring_function: Which docking scoring function to use.
    :param exhaustiveness: Docking exhaustiveness parameter.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    docking_settings = {
        "executable": executable,
        "exhaustiveness": exhaustiveness,
        "scoring_function": scoring_function,
    }

    workflow = stjames.BatchDockingWorkflow(
        initial_smiles_list=smiles_list,
        target=protein,
        protein=protein,
        pocket=pocket,
        docking_settings=docking_settings,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "batch_docking",
        "workflow_data": workflow.model_dump(serialize_as_any=True, mode="json"),
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


__all__ = ["BatchDockingResult", "submit_batch_docking_workflow"]
