"""Pocket detection workflow - detect potential binding sites on a protein."""

from dataclasses import dataclass

import stjames

from rowan.folder import Folder
from rowan.protein import Protein
from rowan.types import ProteinUUID
from rowan.utils import api_client

from .base import Workflow, WorkflowResult, register_result


@dataclass(frozen=True, slots=True)
class Pocket:
    """A detected binding pocket.

    Attributes:
        sphere_centers: centers of detected spheres, in Å
        sphere_radii: radii of detected spheres, in Å
        volume: pocket volume, in Å³
        score: druggability/quality score; larger is better
        pocket_center: center of axis-aligned bounding box, in Å
        pocket_sides: side lengths of axis-aligned bounding box, in Å
        residue_numbers: residue numbers lining the pocket
    """

    sphere_centers: tuple[tuple[float, float, float], ...]
    sphere_radii: tuple[float, ...]
    volume: float
    score: float
    pocket_center: tuple[float, float, float]
    pocket_sides: tuple[float, float, float]
    residue_numbers: tuple[int, ...]


@register_result("pocket_detection")
class PocketDetectionResult(WorkflowResult):
    """Result from a pocket-detection workflow."""

    _stjames_class = stjames.PocketDetectionWorkflow

    def __repr__(self) -> str:
        return f"<PocketDetectionResult pockets={len(self.pockets)}>"

    @property
    def pockets(self) -> list[Pocket]:
        """Detected pockets, in the order returned by the backend."""
        raw = getattr(self._workflow, "pockets", []) or []
        return [
            Pocket(
                sphere_centers=tuple(tuple(c) for c in p.sphere_centers),
                sphere_radii=tuple(p.sphere_radii),
                volume=p.volume,
                score=p.score,
                pocket_center=tuple(p.pocket_center),
                pocket_sides=tuple(p.pocket_sides),
                residue_numbers=tuple(p.residue_numbers),
            )
            for p in raw
        ]


def submit_pocket_detection_workflow(
    protein: Protein | ProteinUUID,
    merge_distance: float = 1.75,
    name: str = "Pocket Detection Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
    is_draft: bool = False,
) -> Workflow[PocketDetectionResult]:
    """Submits a pocket-detection workflow to the API.

    Args:
        protein: protein to analyze. Can be a UUID or a Protein object
        merge_distance: distance for merging pocket spheres, in Å
        name: name of the workflow
        folder_uuid: UUID of the folder to place the workflow in
        folder: destination folder
        max_credits: maximum credits for the workflow
        webhook_url: URL that Rowan will POST to when the workflow completes
        is_draft: save as a draft without starting execution

    Returns:
        workflow object representing the submitted workflow

    Raises:
        httpx.HTTPStatusError: request to the API fails
    """
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder:
        folder_uuid = folder.uuid
    if isinstance(protein, Protein):
        protein = protein.uuid

    workflow = stjames.PocketDetectionWorkflow(
        protein=protein,
        merge_distance=merge_distance,
    )

    data = {
        "workflow_type": "pocket_detection",
        "workflow_data": workflow.model_dump(mode="json"),
        "name": name,
        "folder_uuid": folder_uuid,
        "max_credits": max_credits,
        "webhook_url": webhook_url,
        "is_draft": is_draft,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow[PocketDetectionResult](**response.json())
