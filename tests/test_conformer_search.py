"""Tests for conformer-search workflow submission."""

from collections.abc import Iterator
from contextlib import contextmanager
from unittest.mock import MagicMock

import pytest
import stjames

import rowan
from rowan.workflows import conformer_search as conformer_search_module


def _workflow_response() -> dict[str, object]:
    """Build a minimal workflow API response."""
    return {
        "name": "test",
        "uuid": "workflow-uuid",
        "created_at": "2026-09-15T00:00:00Z",
        "object_status": 0,
        "parent_uuid": "folder-uuid",
        "notes": "",
        "starred": False,
        "public": False,
        "object_type": "conformer_search",
        "object_data": {},
        "email_when_complete": False,
        "credits_charged": 0,
        "object_logfile": "",
    }


def _mock_submission(monkeypatch: pytest.MonkeyPatch) -> MagicMock:
    """Replace the conformer-search API client and return its mock client."""
    client = MagicMock()
    client.post.return_value.json.return_value = _workflow_response()

    @contextmanager
    def mock_api_client() -> Iterator[MagicMock]:
        yield client

    monkeypatch.setattr(conformer_search_module, "api_client", mock_api_client)
    return client


def test_conformer_constraints_are_released_for_ts_refinement(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Keep generation constrained while releasing TS refinement constraints."""
    client = _mock_submission(monkeypatch)
    constraint = rowan.Constraint(constraint_type="freeze_atoms", atoms=[1, 2])
    conf_gen_settings = rowan.iMTDGCSettings(constraints=[constraint])
    optimization_stage = stjames.Settings(
        method=stjames.Method.OMOL25_CONSERVING_S,
        tasks=[stjames.Task.OPTIMIZE],
        opt_settings={
            "transition_state": True,
            "constraints": [constraint],
        },
    )
    multistage_opt_settings = stjames.MultiStageOptSettings(
        optimization_settings=[optimization_stage]
    )

    rowan.submit_conformer_search_workflow(
        initial_molecule=rowan.Molecule.from_xyz("2\n\nH 0 0 0\nH 0 0 0.74"),
        conf_gen_settings=conf_gen_settings,
        multistage_opt_settings=multistage_opt_settings,
    )

    workflow_data = client.post.call_args.kwargs["json"]["workflow_data"]
    assert workflow_data["conf_gen_settings"]["constraints"] == [
        {"constraint_type": "freeze_atoms", "atoms": [1, 2]}
    ]
    assert (
        workflow_data["multistage_opt_settings"]["optimization_settings"][0]["opt_settings"][
            "constraints"
        ]
        == []
    )
    assert optimization_stage.opt_settings.constraints == [constraint]
