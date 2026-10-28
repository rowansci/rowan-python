"""Pocket detection workflow - detect potential binding sites on a protein."""

from dataclasses import dataclass

import stjames
from stjames import MixMDSettings, PocketeerSettings

from rowan.folder import Folder
from rowan.protein import Protein, retrieve_protein
from rowan.types import ProteinUUID
from rowan.utils import api_client

from .base import Message, Workflow, WorkflowResult, parse_messages, register_result


@dataclass(frozen=True, slots=True)
class MixMDPocketDiagnostics:
    """Probe enrichment for a pocket found by mixed-solvent MD.

    Attributes:
        density_center: density-weighted center of the enriched region, in Å
        peak_enrichment: maximum probe density relative to bulk solvent
        integrated_excess_density: probe molecules above bulk, integrated over the region
    """

    density_center: tuple[float, float, float]
    peak_enrichment: float
    integrated_excess_density: float


@dataclass(frozen=True, slots=True)
class Pocket:
    """A detected binding pocket.

    Mixed-solvent MD pockets describe regions of enriched probe density, so their spheres
    and volume outline that region rather than the cavity.

    Attributes:
        sphere_centers: centers of detected spheres, in Å
        sphere_radii: radii of detected spheres, in Å
        volume: pocket volume, in Å³
        score: method-specific ranking score; larger is better
        pocket_center: center of axis-aligned bounding box, in Å
        pocket_sides: side lengths of axis-aligned bounding box, in Å
        residue_numbers: residue numbers lining the pocket
        mixmd: probe-enrichment details for mixed-solvent MD pockets
    """

    sphere_centers: tuple[tuple[float, float, float], ...]
    sphere_radii: tuple[float, ...]
    volume: float
    score: float
    pocket_center: tuple[float, float, float]
    pocket_sides: tuple[float, float, float]
    residue_numbers: tuple[int, ...]
    mixmd: MixMDPocketDiagnostics | None = None


@dataclass(frozen=True, slots=True)
class MixMDReplicaDiagnostics:
    """One mixed-solvent MD replica.

    Attributes:
        index: zero-based replica index
        seed: replica random seed
        num_frames: analyzed frames; zero for a failed replica
        num_probes: probe molecules in the simulation
        largest_probe_cluster_fraction: mean fraction of probes in the largest cluster
        trajectory_uuid: UUID of the stored trajectory
        topology_uuid: UUID of the trajectory topology protein
        density_map_uuid: UUID of the stored probe-density map
    """

    index: int
    seed: int
    num_frames: int
    num_probes: int
    largest_probe_cluster_fraction: float | None
    trajectory_uuid: str | None
    topology_uuid: str | None
    density_map_uuid: str | None


@dataclass(frozen=True, slots=True)
class MixMDResults:
    """Mixed-solvent MD replicas and consensus artifacts.

    Attributes:
        reference_uuid: UUID of the protein defining pocket and map coordinates
        density_map_uuid: UUID of the consensus probe-density map
        replicas: per-replica results, including failed replicas
    """

    reference_uuid: str | None
    density_map_uuid: str | None
    replicas: tuple[MixMDReplicaDiagnostics, ...]


def _mixmd_diagnostics(
    raw: stjames.MixMDPocketDiagnostics | None,
) -> MixMDPocketDiagnostics | None:
    """Convert stjames probe-enrichment diagnostics, if present."""
    if raw is None:
        return None
    return MixMDPocketDiagnostics(
        density_center=tuple(raw.density_center),
        peak_enrichment=raw.peak_enrichment,
        integrated_excess_density=raw.integrated_excess_density,
    )


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
                mixmd=_mixmd_diagnostics(p.mixmd),
            )
            for p in raw
        ]

    @property
    def mixmd(self) -> MixMDResults | None:
        """Mixed-solvent MD replicas and artifacts, or `None` for Pocketeer runs."""
        if not (raw := getattr(self._workflow, "mixmd", None)):
            return None
        return MixMDResults(
            reference_uuid=raw.reference_uuid,
            density_map_uuid=raw.density_map_uuid,
            replicas=tuple(
                MixMDReplicaDiagnostics(
                    index=r.index,
                    seed=r.seed,
                    num_frames=r.num_frames,
                    num_probes=r.num_probes,
                    largest_probe_cluster_fraction=r.largest_probe_cluster_fraction,
                    trajectory_uuid=r.trajectory_uuid,
                    topology_uuid=r.topology_uuid,
                    density_map_uuid=r.density_map_uuid,
                )
                for r in raw.replicas
            ),
        )

    @property
    def messages(self) -> list[Message]:
        """Any messages or warnings from the workflow, including MixMD reliability warnings."""
        return parse_messages(getattr(self._workflow, "messages", None))

    def get_reference_protein(self) -> Protein | None:
        """Fetch and cache the protein defining mixed-solvent MD pocket and map coordinates.

        Returns:
            protein object, or `None` for Pocketeer runs
        """
        if not (mixmd := self.mixmd) or not mixmd.reference_uuid:
            return None
        if "reference_protein" not in self._cache:
            self._cache["reference_protein"] = retrieve_protein(
                mixmd.reference_uuid, workflow_uuid=self.workflow_uuid
            )
        return self._cache["reference_protein"]


def submit_pocket_detection_workflow(
    protein: Protein | ProteinUUID,
    merge_distance: float | None = None,
    settings: PocketeerSettings | MixMDSettings | None = None,
    validate_forcefield: bool = True,
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
        merge_distance: Pocketeer distance for merging pocket spheres, in Å; defaults to 1.75.
            Set it on `PocketeerSettings` instead when passing `settings`
        settings: method settings; `PocketeerSettings` (default) for geometric detection or
            `MixMDSettings` for mixed-solvent molecular dynamics
        validate_forcefield: validate the protein forcefield compatibility before submitting a
            mixed-solvent MD run. Raises an error early if the protein cannot be parameterized
            or has clashing residues
        name: name of the workflow
        folder_uuid: UUID of the folder to place the workflow in
        folder: destination folder
        max_credits: maximum credits for the workflow
        webhook_url: URL that Rowan will POST to when the workflow completes
        is_draft: save as a draft without starting execution

    Returns:
        workflow object representing the submitted workflow

    Raises:
        ValueError: both `merge_distance` and `settings`, or both `folder` and `folder_uuid`,
            are given
        httpx.HTTPStatusError: request to the API fails
    """
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if settings is not None and merge_distance is not None:
        raise ValueError("Set `merge_distance` on `PocketeerSettings` when passing `settings`.")
    if folder:
        folder_uuid = folder.uuid
    if isinstance(protein, Protein):
        protein = protein.uuid

    if settings is None:
        settings = (
            PocketeerSettings()
            if merge_distance is None
            else PocketeerSettings(merge_distance=merge_distance)
        )
    if validate_forcefield and isinstance(settings, MixMDSettings):
        Protein(uuid=protein).validate_protein_forcefield()

    workflow = stjames.PocketDetectionWorkflow(protein=protein, settings=settings)

    data = {
        "workflow_type": "pocket_detection",
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
        return Workflow[PocketDetectionResult](**response.json())
