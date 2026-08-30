"""Tests for the docking workflow, focused on induced-fit docking support."""

from collections.abc import Iterator
from contextlib import contextmanager
from typing import Any
from unittest.mock import MagicMock

import stjames
from pytest import MonkeyPatch, approx, raises

import rowan
from rowan.workflows import docking as docking_module


def _docking_workflow_data(*, induced_fit: bool = False, **score_kwargs: Any) -> dict[str, Any]:
    """Build a minimal, schema-valid DockingWorkflow response with one score."""
    wf = stjames.DockingWorkflow(
        initial_molecule=stjames.Molecule(
            atoms=[stjames.Atom(atomic_number=8, position=(0.0, 0.0, 0.0))],
            charge=0,
            multiplicity=1,
        ),
        protein="11111111-1111-1111-1111-111111111111",
        pocket=((0.0, 0.0, 0.0), (10.0, 10.0, 10.0)),
        docking_settings=stjames.VinaSettings(),
        induced_fit_settings=stjames.InducedFitSettings() if induced_fit else None,
        scores=[
            stjames.Score(
                pose="33333333-3333-3333-3333-333333333333",
                complex_pdb=None,
                score=-7.5,
                posebusters_valid=True,
                **score_kwargs,
            )
        ],
    )
    return wf.model_dump(serialize_as_any=True, mode="json")


def _docking_result(
    *, induced_fit: bool = False, **score_kwargs: Any
) -> docking_module.DockingResult:
    """Build a parsed docking result with one score."""
    data = _docking_workflow_data(induced_fit=induced_fit, **score_kwargs)
    return docking_module.DockingResult(
        workflow_data=data,
        workflow_type="docking",
        workflow_uuid="workflow-uuid",
    )


def test_submit_docking_workflow_posts_induced_fit_settings(monkeypatch: MonkeyPatch) -> None:
    """Induced-fit settings reach the workflow payload on DockingWorkflow, not VinaSettings."""
    client = MagicMock()
    client.post.return_value.json.return_value = {
        "name": "docking",
        "uuid": "workflow-uuid",
        "created_at": "2026-07-17T00:00:00Z",
        "object_status": 0,
        "parent_uuid": "folder-uuid",
        "notes": "",
        "starred": False,
        "public": False,
        "object_type": "docking",
        "object_data": {},
        "email_when_complete": False,
        "credits_charged": 0,
        "object_logfile": "",
    }

    @contextmanager
    def mock_api_client() -> Iterator[MagicMock]:
        yield client

    monkeypatch.setattr(docking_module, "api_client", mock_api_client)

    docking_module.submit_docking_workflow(
        protein="protein-uuid",
        pocket=[[0.0, 0.0, 0.0], [10.0, 10.0, 10.0]],
        initial_molecule=stjames.Molecule(
            atoms=[stjames.Atom(atomic_number=8, position=(0.0, 0.0, 0.0))],
            charge=0,
            multiplicity=1,
        ),
        induced_fit_settings=rowan.InducedFitSettings(max_receptors=3),
    )

    posted = client.post.call_args.kwargs["json"]
    assert posted["workflow_data"]["induced_fit_settings"] == {
        "max_receptors": 3,
        "flexible_sidechain_radius": 5.0,
    }
    assert "induced_fit_settings" not in posted["workflow_data"]["docking_settings"]
    assert rowan.InducedFitSettings is stjames.InducedFitSettings


def test_docking_result_scores_expose_induced_fit_fields() -> None:
    """The composite score, receptor strain, and geometry penalty round-trip onto DockingScore."""
    result = _docking_result(
        induced_fit=True,
        receptor_strain=1.234,
        geometry_penalty=0.0,
        induced_fit_score=-6.266,
        induced_receptor_pdb="22222222-2222-2222-2222-222222222222",
    )

    [score] = result.scores
    assert score.score == approx(-7.5)
    assert score.receptor_strain == approx(1.234)
    assert score.geometry_penalty == approx(0.0)
    assert score.induced_fit_score == approx(-6.266)
    assert score.induced_receptor_pdb == "22222222-2222-2222-2222-222222222222"
    assert repr(result) == "<DockingResult poses=1 best_score=-6.266>"


def test_docking_result_scores_defaults_induced_fit_fields_to_none() -> None:
    """A rigid pose without induced-fit fields parses with them all unset."""
    result = _docking_result()

    [score] = result.scores
    assert score.receptor_strain is None
    assert score.geometry_penalty is None
    assert score.induced_fit_score is None
    assert score.induced_receptor_pdb is None


def test_get_induced_receptor_fetches_by_uuid(monkeypatch: MonkeyPatch) -> None:
    """The induced receptor is fetched via its own UUID, cached, and reused."""
    result = _docking_result(
        induced_fit=True,
        induced_receptor_pdb="22222222-2222-2222-2222-222222222222",
    )
    retrieve = MagicMock(return_value="induced-receptor")
    monkeypatch.setattr(docking_module, "retrieve_protein", retrieve)

    receptor = result.get_induced_receptor(0)
    result.get_induced_receptor(0)

    assert receptor == "induced-receptor"
    retrieve.assert_called_once_with(
        "22222222-2222-2222-2222-222222222222", workflow_uuid="workflow-uuid"
    )


def test_get_induced_receptor_raises_without_uuid() -> None:
    """A rigid pose has no induced receptor to fetch."""
    result = _docking_result()

    with raises(ValueError, match="no induced receptor"):
        result.get_induced_receptor(0)
