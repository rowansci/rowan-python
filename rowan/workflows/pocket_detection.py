"""Pocket detection workflow - detect potential binding sites on a protein."""

from dataclasses import dataclass

import stjames

from ..folder import Folder
from ..protein import Protein
from ..utils import api_client
from .base import Workflow, WorkflowResult, register_result


@dataclass(frozen=True, slots=True)
class Pocket:
    """
    A detected binding pocket.

    :param sphere_centers: centers of detected spheres, in Å
    :param sphere_radii: radii of detected spheres, in Å
    :param volume: pocket volume, in Å³
    :param score: druggability/quality score; larger is better
    :param pocket_center: center of axis-aligned bounding box, in Å
    :param pocket_sides: side lengths of axis-aligned bounding box, in Å
    :param residue_numbers: residue numbers lining the pocket
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
    protein: str | Protein,
    merge_distance: float = 1.75,
    name: str = "Pocket Detection Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
    is_draft: bool = False,
) -> Workflow:
    """
    Submits a pocket-detection workflow to the API.

    :param protein: protein to analyze. Can be a UUID or a Protein object.
    :param merge_distance: distance for merging pocket spheres, in Å
    :param name: name of the workflow
    :param folder_uuid: UUID of the folder to place the workflow in
    :param folder: Folder object to store the workflow in
    :param max_credits: maximum number of credits to use for the workflow
    :param webhook_url: URL that Rowan will POST to when the workflow completes
    :param is_draft: if True, submit the workflow as a draft without starting execution
    :returns: Workflow object representing the submitted workflow
    :raises requests.HTTPError: if the request to the API fails
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
        return Workflow(**response.json())
