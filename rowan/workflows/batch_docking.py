"""Batch docking workflow - high-throughput molecular docking."""

import stjames

from rowan.folder import Folder
from rowan.protein import Protein
from rowan.types import ProteinUUID
from rowan.utils import api_client

from .base import Workflow, WorkflowResult, register_result
from .docking import DockingScore


@register_result("batch_docking")
class BatchDockingResult(WorkflowResult):
    """Result from a batch-docking workflow."""

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
        scores = self._workflow.best_scores or []
        padded = list(scores) + [None] * (len(smiles_list) - len(scores))
        return dict(zip(smiles_list, padded, strict=True))

    @property
    def refined_scores(self) -> dict[str, DockingScore | None]:
        """Saved docking results indexed by input SMILES."""
        smiles_list = self._workflow.initial_smiles_list
        scores = self._workflow.refined_scores or []
        padded = list(scores) + [None] * (len(smiles_list) - len(scores))
        return {
            smiles: (
                DockingScore(
                    score=score.score,
                    pose=score.pose,
                    complex_pdb=score.complex_pdb,
                    posebusters_valid=score.posebusters_valid,
                    strain=score.strain,
                    rmsd=score.rmsd,
                    mmgbsa_score=score.mmgbsa_score,
                )
                if score is not None
                else None
            )
            for smiles, score in zip(smiles_list, padded, strict=True)
        }


def submit_batch_docking_workflow(
    smiles_list: list[str],
    protein: Protein | ProteinUUID,
    pocket: list[list[float]],
    executable: str = "vina",
    scoring_function: str = "vinardo",
    exhaustiveness: float = 8,
    num_poses_to_save: int = 0,
    run_mmgbsa: bool = False,
    name: str = "Batch Docking Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
    is_draft: bool = False,
) -> Workflow[BatchDockingResult]:
    """Submits a batch-docking workflow to the API.

    Args:
        smiles_list: SMILES strings to dock
        protein: protein to dock (UUID or Protein object)
        pocket: binding pocket as `[[cx, cy, cz], [sx, sy, sz]]` – center (Å) and box size (Å)
        executable: which docking implementation to use
        scoring_function: which docking scoring function to use
        exhaustiveness: docking exhaustiveness parameter
        num_poses_to_save: number of top-scoring compounds whose best pose to save
        run_mmgbsa: whether to refine the saved poses with MM/GBSA. Ignored when
            `num_poses_to_save` is zero
        name: name of the workflow
        folder_uuid: UUID of the folder to place the workflow in
        folder: destination folder
        max_credits: maximum number of credits to use
        webhook_url: URL that Rowan will POST to when the workflow completes
        is_draft: save as a draft without starting execution

    Returns:
        submitted workflow

    Raises:
        httpx.HTTPStatusError: request to the API fails
    """
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder:
        folder_uuid = folder.uuid
    if isinstance(protein, Protein):
        protein = protein.uuid
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
        num_poses_to_save=num_poses_to_save,
        run_mmgbsa=run_mmgbsa,
    )

    data = {
        "workflow_type": "batch_docking",
        "workflow_data": workflow.model_dump(serialize_as_any=True, mode="json"),
        "name": name,
        "folder_uuid": folder_uuid,
        "max_credits": max_credits,
        "webhook_url": webhook_url,
        "is_draft": is_draft,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow[BatchDockingResult](**response.json())
