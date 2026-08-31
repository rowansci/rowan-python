"""Regression tests for docking workflow updates."""

import inspect
from collections.abc import Iterator
from contextlib import contextmanager
from unittest.mock import MagicMock

import stjames
from pytest import MonkeyPatch

from rowan.workflows import analogue_docking as analogue_docking_module
from rowan.workflows import batch_docking as batch_docking_module
from rowan.workflows.batch_docking import BatchDockingResult


def _workflow_response() -> dict[str, object]:
    """Build a minimal workflow API response."""
    return {
        "name": "test",
        "uuid": "workflow-uuid",
        "created_at": "2026-08-24T00:00:00Z",
        "object_status": 0,
        "parent_uuid": "folder-uuid",
        "notes": "",
        "starred": False,
        "public": False,
        "object_type": "batch_docking",
        "object_data": {},
        "email_when_complete": False,
        "credits_charged": 0,
        "object_logfile": "",
    }


def _mock_submission(monkeypatch: MonkeyPatch) -> MagicMock:
    """Replace the batch docking API client and return its mock client."""
    client = MagicMock()
    client.post.return_value.json.return_value = _workflow_response()

    @contextmanager
    def mock_api_client() -> Iterator[MagicMock]:
        yield client

    monkeypatch.setattr(batch_docking_module, "api_client", mock_api_client)
    return client


def test_batch_docking_submission_controls(monkeypatch: MonkeyPatch) -> None:
    """Submit default and opt-in pose-saving and refinement controls."""
    client = _mock_submission(monkeypatch)

    batch_docking_module.submit_batch_docking_workflow(
        ["CC", "CCC"],
        "protein-uuid",
        [[0.0, 0.0, 0.0], [10.0, 10.0, 10.0]],
    )
    payload = client.post.call_args.kwargs["json"]["workflow_data"]
    assert payload["num_poses_to_save"] == 0
    assert payload["run_mmgbsa"] is False

    batch_docking_module.submit_batch_docking_workflow(
        ["CC", "CCC"],
        "protein-uuid",
        [[0.0, 0.0, 0.0], [10.0, 10.0, 10.0]],
        num_poses_to_save=2,
        run_mmgbsa=True,
    )
    payload = client.post.call_args.kwargs["json"]["workflow_data"]
    assert payload["num_poses_to_save"] == 2
    assert payload["run_mmgbsa"] is True


def test_batch_docking_results_align_with_input_smiles() -> None:
    """Align partial raw and refined scores with their input SMILES."""

    workflow = stjames.BatchDockingWorkflow(
        protein="protein-uuid",
        target="protein-uuid",
        initial_smiles_list=["CC", "CCC", "CCCC"],
        pocket=((0, 0, 0), (10, 10, 10)),
        best_scores=[-5.0, None],
        refined_scores=[
            {
                "pose": "pose-uuid",
                "complex_pdb": "complex-uuid",
                "score": -5.0,
                "posebusters_valid": True,
                "mmgbsa_score": -20.0,
            },
            None,
        ],
    )
    workflow_data = workflow.model_dump(mode="json")
    result = BatchDockingResult(
        workflow_data=workflow_data,
        workflow_type="batch_docking",
        workflow_uuid="workflow-uuid",
    )
    assert result.scores == {"CC": -5.0, "CCC": None, "CCCC": None}
    refined = result.refined_scores["CC"]
    assert refined is not None
    assert refined.pose == "pose-uuid"
    assert refined.complex_pdb == "complex-uuid"
    assert refined.mmgbsa_score == -20.0
    assert result.refined_scores["CCC"] is None
    assert result.refined_scores["CCCC"] is None


def test_analogue_docking_defaults_match_stjames() -> None:
    """Keep public analogue controls aligned with StJames."""
    parameters = inspect.signature(
        analogue_docking_module.submit_analogue_docking_workflow
    ).parameters
    assert parameters["num_conformers_per_analogue"].default == 20
    assert "require_posebusters" not in parameters
    assert stjames.AnalogueDockingWorkflow.model_fields["num_conformers_per_analogue"].default == 20
    assert stjames.AnalogueDockingWorkflow.model_fields["require_posebusters"].default is False
