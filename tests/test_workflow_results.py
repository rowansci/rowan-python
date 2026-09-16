"""Tests for result retrieval and streaming from typed workflows."""

from unittest.mock import MagicMock

import stjames
from pytest import MonkeyPatch, fixture, raises

import rowan
from rowan.workflows import admet, base


@fixture
def workflow_client(monkeypatch: MonkeyPatch) -> MagicMock:
    """Mock an ADMET workflow submission and completed response."""
    client = MagicMock()
    client.__enter__.return_value = client
    response = {
        "name": "ethanol",
        "uuid": "workflow-uuid",
        "created_at": "2026-09-15T00:00:00Z",
        "object_status": stjames.Status.COMPLETED_OK,
        "parent_uuid": "folder-uuid",
        "notes": "",
        "starred": False,
        "public": False,
        "object_type": "admet",
        "object_data": {"initial_smiles": "CCO", "properties": {"molecular_weight": 46.07}},
        "email_when_complete": False,
        "credits_charged": 0,
        "object_logfile": "",
    }
    client.post.return_value.json.return_value = response | {"object_status": stjames.Status.QUEUED}
    client.get.return_value.json.return_value = response
    monkeypatch.setattr(base, "api_client", lambda: client)
    monkeypatch.setattr(admet, "api_client", lambda: client)
    return client


def test_submitted_workflow_refresh_and_stream(workflow_client: MagicMock) -> None:
    """Refresh a submitted workflow and stream its completed properties."""
    workflow = rowan.submit_admet_workflow("CCO")
    refreshed = workflow.fetch_latest()
    assert workflow.status == stjames.Status.QUEUED
    assert refreshed.status == stjames.Status.COMPLETED_OK
    results = list(refreshed.stream_result(poll_interval=0))
    assert len(results) == 1
    assert results[0].properties == {"molecular_weight": 46.07}
    assert refreshed.result().properties == {"molecular_weight": 46.07}


def test_retrieve_expected_result(workflow_client: MagicMock) -> None:
    """Retrieve and parse the expected workflow result."""
    workflow = rowan.retrieve_workflow("workflow-uuid", result_type=rowan.ADMETResult)
    assert workflow.result().properties == {"molecular_weight": 46.07}


def test_retrieve_wrong_result(workflow_client: MagicMock) -> None:
    """Reject a UUID belonging to another workflow kind."""
    with raises(ValueError):
        rowan.retrieve_workflow("workflow-uuid", result_type=rowan.DockingResult)


def test_retrieve_without_result_type(workflow_client: MagicMock) -> None:
    """Keep untyped retrieval compatible with registered results."""
    workflow = rowan.retrieve_workflow("workflow-uuid")
    assert workflow.result().data["properties"] == {"molecular_weight": 46.07}
