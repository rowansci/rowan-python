"""Tests for tautomer-search workflow submission."""

from collections.abc import Iterator
from contextlib import contextmanager
from unittest.mock import MagicMock

from pytest import MonkeyPatch

import rowan
from rowan.workflows import tautomer_search as tautomer_module


def _mock_submission(monkeypatch: MonkeyPatch) -> MagicMock:
    """Replace the workflow API client and return its mock client."""
    client = MagicMock()
    client.post.return_value.json.return_value = {
        "name": "test",
        "uuid": "workflow-uuid",
        "created_at": "2026-09-01T00:00:00Z",
        "object_status": 0,
        "parent_uuid": "folder-uuid",
        "notes": "",
        "starred": False,
        "public": False,
        "object_type": "tautomers",
        "object_data": {},
        "email_when_complete": False,
        "credits_charged": 0,
        "object_logfile": "",
    }

    @contextmanager
    def mock_api_client() -> Iterator[MagicMock]:
        yield client

    monkeypatch.setattr(tautomer_module, "api_client", mock_api_client)
    return client


def test_submit_tautomer_search_screening_window(monkeypatch: MonkeyPatch) -> None:
    """Serialize the default and overridden tautomer screening windows."""
    client = _mock_submission(monkeypatch)
    molecule = rowan.Molecule.from_xyz("2\nhydrogen\nH 0 0 0\nH 0 0 0.74")

    tautomer_module.submit_tautomer_search_workflow(molecule)
    workflow_data = client.post.call_args.kwargs["json"]["workflow_data"]
    assert workflow_data["screening_window"] == 10.0

    tautomer_module.submit_tautomer_search_workflow(molecule, screening_window=25.0)
    workflow_data = client.post.call_args.kwargs["json"]["workflow_data"]
    assert workflow_data["screening_window"] == 25.0
