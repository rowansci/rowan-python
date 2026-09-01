"""Tests for workflow submission requests."""

from contextlib import contextmanager
from typing import Any
from unittest.mock import MagicMock, call

import pytest
import stjames
from pytest import MonkeyPatch

from rowan.workflows import base
from rowan.workflows import solubility as solubility_module


def _workflow_response(index: int) -> dict[str, Any]:
    """Build one minimal workflow API response."""
    return {
        "name": f"workflow-{index}",
        "uuid": f"workflow-uuid-{index}",
        "created_at": "2026-07-17T00:00:00Z",
        "object_status": 0,
        "parent_uuid": "folder-uuid",
        "notes": "",
        "starred": False,
        "public": False,
        "object_type": "solubility",
        "object_data": {},
        "email_when_complete": False,
        "credits_charged": 0,
        "object_logfile": "",
        "submission_group_uuid": "submission-group-uuid",
    }


def test_submit_solubility_workflow_group_posts_one_group(monkeypatch: MonkeyPatch) -> None:
    """Submit validated solubility settings in one group request."""
    client = MagicMock()
    client.post.return_value.json.return_value = [_workflow_response(1), _workflow_response(2)]

    @contextmanager
    def mock_api_client():
        yield client

    monkeypatch.setattr(base, "api_client", mock_api_client)

    workflows = solubility_module.submit_solubility_workflow_group(
        ["CCO", "CCN"],
        method="kingfisher",
        names=["ethanol", "ethylamine"],
        folder_uuid="folder-uuid",
    )

    client.post.assert_called_once_with(
        "/workflow/submit_group",
        json={
            "workflow_type": "solubility",
            "workflow_data": {
                "solubility_method": "kingfisher",
                "solvents": ["O"],
                "temperatures": [298.15],
            },
            "initial_smileses": ["CCO", "CCN"],
            "names": ["ethanol", "ethylamine"],
            "folder_uuid": "folder-uuid",
            "max_credits": None,
            "webhook_url": None,
        },
    )
    assert len(workflows) == 2
    assert {workflow.submission_group_uuid for workflow in workflows} == {"submission-group-uuid"}


def test_batch_submit_workflow_preserves_individual_submission_behavior(
    monkeypatch: MonkeyPatch,
) -> None:
    """Keep the established batch helper as repeated individual submissions."""
    submit = MagicMock(side_effect=["first", "second"])
    monkeypatch.setattr(base, "submit_workflow", submit)

    result = base.batch_submit_workflow(
        workflow_type="solubility",
        workflow_data={"solubility_method": "fastsolv"},
        initial_smileses=["CCO", "CCN"],
        names=["ethanol", "ethylamine"],
        folder_uuid="folder-uuid",
    )

    assert result == ["first", "second"]
    assert submit.call_count == 2


def test_retrieve_workflows_uses_bounded_batches(monkeypatch: MonkeyPatch) -> None:
    """Split large full-workflow retrievals into bounded API requests."""
    client = MagicMock()
    first_response = MagicMock()
    first_response.json.return_value = [_workflow_response(index) for index in range(100)]
    second_response = MagicMock()
    second_response.json.return_value = [_workflow_response(100)]
    client.post.side_effect = [first_response, second_response]

    @contextmanager
    def mock_api_client():
        yield client

    monkeypatch.setattr(base, "api_client", mock_api_client)

    uuids = [f"workflow-uuid-{index}" for index in range(101)]
    workflows = base.retrieve_workflows(uuids)

    assert len(workflows) == 101
    assert client.post.call_args_list == [
        call(
            "/workflow/batch_retrieve",
            json={"uuids": uuids[:100]},
        ),
        call(
            "/workflow/batch_retrieve",
            json={"uuids": uuids[100:]},
        ),
    ]


@pytest.mark.parametrize("status", [stjames.Status.FAILED, stjames.Status.STOPPED])
def test_workflow_result_includes_log_for_unsuccessful_workflow(
    monkeypatch: MonkeyPatch,
    status: stjames.Status,
) -> None:
    """Include the refreshed workflow log in failed and stopped result errors."""
    logfile = "Calculation failed because the input was invalid."

    def add_logfile(workflow: base.Workflow, in_place: bool = True) -> base.Workflow:
        assert in_place
        workflow.logfile = logfile
        return workflow

    monkeypatch.setattr(base.Workflow, "fetch_latest", add_logfile)
    response = _workflow_response(1) | {"object_status": status}
    workflow = base.Workflow.model_validate(response)

    with pytest.raises(base.WorkflowError) as error_info:
        workflow.result(wait=False)

    assert error_info.value.logfile == logfile
    assert str(error_info.value) == (
        f"Workflow 'workflow-1' {status.name.lower()} (uuid=workflow-uuid-1). "
        "See WorkflowError.logfile for diagnostic details."
    )
