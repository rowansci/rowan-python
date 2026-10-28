"""Tests for pocket detection submission and mixed-solvent MD results."""

from collections.abc import Iterator
from contextlib import contextmanager
from typing import Any
from unittest.mock import MagicMock

import pytest
from stjames.workflows import pocket_detection as stjames_pocket_detection

import rowan
from rowan.workflows import pocket_detection
from rowan.workflows.pocket_detection import MixMDPocketDiagnostics, PocketDetectionResult

REFERENCE_UUID = "3f2b8c1e-5d4a-4b6e-9f10-2a3b4c5d6e7f"
MAP_UUID = "8a7b6c5d-4e3f-4a1b-9c8d-7e6f5a4b3c2d"


@pytest.fixture
def client(monkeypatch: pytest.MonkeyPatch) -> MagicMock:
    """Replace the pocket-detection API client and return the mock."""
    client = MagicMock()
    client.post.return_value.json.return_value = {
        "name": "test",
        "uuid": "workflow-uuid",
        "created_at": "2026-10-22T00:00:00Z",
        "object_status": 0,
        "parent_uuid": "folder-uuid",
        "notes": "",
        "starred": False,
        "public": False,
        "object_type": "pocket_detection",
        "object_data": {},
        "email_when_complete": False,
        "credits_charged": 0,
        "object_logfile": "",
    }

    @contextmanager
    def mock_api_client() -> Iterator[MagicMock]:
        yield client

    monkeypatch.setattr(pocket_detection, "api_client", mock_api_client)
    return client


def _submitted_settings(client: MagicMock) -> dict[str, Any]:
    """Return the settings sent in the workflow submission."""
    return client.post.call_args.kwargs["json"]["workflow_data"]["settings"]


@pytest.mark.parametrize(("merge_distance", "expected"), [(None, 1.75), (2.5, 2.5)])
def test_pocketeer_submission_sends_merge_distance_in_settings(
    client: MagicMock, merge_distance: float | None, expected: float
) -> None:
    """Send Pocketeer by default with the requested merge distance."""
    rowan.submit_pocket_detection_workflow("protein-uuid", merge_distance=merge_distance)
    assert _submitted_settings(client) == {"settings_type": "pocketeer", "merge_distance": expected}


def test_mixmd_submission_validates_forcefield(
    client: MagicMock, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Validate the protein, then send MixMD settings with their defaults."""
    validate = MagicMock()
    monkeypatch.setattr(rowan.Protein, "validate_protein_forcefield", validate)
    settings = rowan.MixMDSettings(probe=rowan.Solvent.ISOPROPANOL, num_trajectories=5)
    rowan.submit_pocket_detection_workflow("protein-uuid", settings=settings)

    validate.assert_called_once()
    submitted = _submitted_settings(client)
    assert submitted["settings_type"] == "mixmd"
    assert submitted["probe"] == rowan.Solvent.ISOPROPANOL.value
    assert submitted["num_trajectories"] == 5
    assert submitted["enrichment_threshold"] == 5


def test_merge_distance_with_settings_is_rejected() -> None:
    """Reject a top-level merge distance alongside explicit settings."""
    with pytest.raises(ValueError, match="merge_distance"):
        rowan.submit_pocket_detection_workflow(
            "protein-uuid", merge_distance=2, settings=rowan.PocketeerSettings()
        )


def test_mixmd_result_exposes_pocket_and_replica_details() -> None:
    """Read probe enrichment, replica diagnostics, and artifact UUIDs, including failures."""
    workflow = stjames_pocket_detection.PocketDetectionWorkflow(
        protein="protein-uuid",
        settings=rowan.MixMDSettings(),
        pockets=[
            stjames_pocket_detection.Pocket(
                sphere_centers=[(1, 2, 3)],
                sphere_radii=[2.2],
                volume=12,
                score=4.5,
                pocket_center=(1, 2, 3),
                pocket_sides=(15, 15, 15),
                residue_numbers=[7, 9],
                mixmd=stjames_pocket_detection.MixMDPocketDiagnostics(
                    density_center=(1.5, 2, 3), peak_enrichment=8, integrated_excess_density=4.5
                ),
            )
        ],
        mixmd=stjames_pocket_detection.MixMDResults(
            reference_uuid=REFERENCE_UUID,
            density_map_uuid=MAP_UUID,
            replicas=[
                stjames_pocket_detection.MixMDReplicaDiagnostics(
                    index=0, seed=11, num_frames=1500, num_probes=56, density_map_uuid=MAP_UUID
                ),
                stjames_pocket_detection.MixMDReplicaDiagnostics(index=1, seed=12),
            ],
        ),
    )
    result = PocketDetectionResult(
        workflow_data=workflow.model_dump(mode="json"),
        workflow_type="pocket_detection",
        workflow_uuid="workflow-uuid",
    )

    assert result.pockets[0].mixmd == MixMDPocketDiagnostics(
        density_center=(1.5, 2, 3), peak_enrichment=8, integrated_excess_density=4.5
    )
    assert result.mixmd is not None
    assert result.mixmd.reference_uuid == REFERENCE_UUID
    assert result.mixmd.replicas[0].density_map_uuid == MAP_UUID
    assert result.mixmd.replicas[1].num_frames == 0
    assert result.mixmd.replicas[1].trajectory_uuid is None


def test_pocketeer_result_has_no_mixmd_details() -> None:
    """Leave MixMD fields empty for geometric pocket detection."""
    workflow = stjames_pocket_detection.PocketDetectionWorkflow(protein="protein-uuid")
    result = PocketDetectionResult(
        workflow_data=workflow.model_dump(mode="json"),
        workflow_type="pocket_detection",
        workflow_uuid="workflow-uuid",
    )
    assert result.mixmd is None
    assert result.get_reference_protein() is None
