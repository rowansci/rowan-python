"""Tests for SepTop graph submission and topology result handling."""

import json
from typing import Literal, TypedDict
from unittest.mock import MagicMock

from pydantic import ValidationError
from pytest import MonkeyPatch, mark, raises

import rowan
from rowan.workflows import rbfe_graph
from stjames import RBFEGraphWorkflow


class _GraphSubmissionOptions(TypedDict, total=False):
    separated_topology_planning_score: float | Literal["auto"] | None
    greedy_scoring: Literal["best", "jaccard", "dummy_atoms"] | None
    mode: Literal["greedy", "star_map"]
    hub_compound_id: str
    generate_intermediate_ligands: bool


@mark.parametrize(
    ("options", "expected", "expected_scoring"),
    [
        ({}, 25, "dummy_atoms"),
        ({"separated_topology_planning_score": "auto"}, 25, "dummy_atoms"),
        ({"greedy_scoring": "dummy_atoms"}, 25, "dummy_atoms"),
        ({"greedy_scoring": "best"}, None, "best"),
        ({"greedy_scoring": "jaccard"}, None, "jaccard"),
        ({"mode": "star_map", "hub_compound_id": "a"}, None, "best"),
        ({"generate_intermediate_ligands": True}, None, "best"),
        (
            {"generate_intermediate_ligands": True, "separated_topology_planning_score": 25},
            25,
            "dummy_atoms",
        ),
        ({"separated_topology_planning_score": None}, None, "best"),
        ({"separated_topology_planning_score": 0}, 0, "dummy_atoms"),
        ({"separated_topology_planning_score": 40.5}, 40.5, "dummy_atoms"),
        ({"separated_topology_planning_score": None, "greedy_scoring": "jaccard"}, None, "jaccard"),
        (
            {"separated_topology_planning_score": None, "mode": "star_map", "hub_compound_id": "a"},
            None,
            "best",
        ),
    ],
)
def test_graph_submission_serializes_planning_score(
    monkeypatch: MonkeyPatch,
    options: _GraphSubmissionOptions,
    expected: float | None,
    expected_scoring: str,
) -> None:
    """Requests preserve planning defaults, opt-out, and compatible scoring modes."""
    client = MagicMock()
    context = MagicMock()
    context.__enter__.return_value = client
    monkeypatch.setattr(rbfe_graph, "api_client", MagicMock(return_value=context))
    monkeypatch.setattr(rbfe_graph, "Workflow", MagicMock())
    molecule = rowan.Molecule.from_xyz("2\nhydrogen\nH 0 0 0\nH 0 0 0.74")
    rowan.submit_relative_binding_free_energy_graph_workflow(
        ligands={"a": molecule, "b": molecule}, **options
    )
    payload = client.post.call_args.kwargs["json"]["workflow_data"]
    assert payload.get("separated_topology_planning_score") == expected
    saved = RBFEGraphWorkflow.model_validate_json(json.dumps(payload))
    restored = RBFEGraphWorkflow.model_validate(saved.model_dump(mode="json"))
    assert restored.separated_topology_planning_score == expected
    assert restored.greedy_scoring == expected_scoring
    assert "septop_threshold" not in payload


@mark.parametrize(
    "options",
    [
        {"greedy_scoring": "best", "separated_topology_planning_score": 25},
        {"greedy_scoring": "jaccard", "separated_topology_planning_score": 25},
        {"mode": "star_map", "hub_compound_id": "a", "separated_topology_planning_score": 25},
        {"separated_topology_planning_score": -1},
    ],
)
def test_graph_submission_rejects_incompatible_planning(
    monkeypatch: MonkeyPatch, options: _GraphSubmissionOptions
) -> None:
    """Invalid planning settings fail before submitting an API request."""
    client = MagicMock()
    monkeypatch.setattr(rbfe_graph, "api_client", client)
    molecule = rowan.Molecule.from_xyz("2\n\nH 0 0 0\nH 0 0 0.74")
    with raises(ValidationError):
        rowan.submit_relative_binding_free_energy_graph_workflow(
            ligands={"a": molecule, "b": molecule}, **options
        )
    client.assert_not_called()


@mark.parametrize(
    ("result_class", "workflow_type"),
    [
        (rowan.RelativeBindingFreeEnergyGraphResult, "rbfe_graph"),
        (
            rowan.RelativeBindingFreeEnergyPerturbationResult,
            "relative_binding_free_energy_perturbation",
        ),
    ],
)
def test_result_preserves_edge_topology(
    result_class: type[rowan.RelativeBindingFreeEnergyGraphResult]
    | type[rowan.RelativeBindingFreeEnergyPerturbationResult],
    workflow_type: str,
) -> None:
    """Graph and FEP result convenience edges preserve the saved SepTop strategy."""
    molecule = rowan.Molecule.from_xyz("2\nhydrogen\nH 0 0 0\nH 0 0 0.74").to_stjames()
    graph = rowan.RBFEGraph(
        edges=[
            rowan.RBFEGraphEdge(
                ligand_a="a",
                ligand_b="b",
                topology=rowan.RBFEEdgeTopology.SEPARATED_TOPOLOGY,
                score=26,
            )
        ]
    )
    result = result_class(
        workflow_data={
            "ligands": {"a": molecule, "b": molecule},
            "graph": graph,
            "protein": "00000000-0000-0000-0000-000000000001",
            "settings": {},
        },
        workflow_type=workflow_type,
        workflow_uuid="graph-uuid",
    )
    assert result.edges[0].topology == rowan.RBFEEdgeTopology.SEPARATED_TOPOLOGY
    assert result.edges[0].core is None
    assert result.graph == graph
