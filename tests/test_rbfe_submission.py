"""Tests for RBFE cofactor submission."""

from collections.abc import Iterator
from contextlib import contextmanager
from unittest.mock import MagicMock, call

import pytest

import rowan
from rowan.workflows import relative_binding_free_energy_perturbation as rbfe
from rowan.workflows.rbfe_graph import RelativeBindingFreeEnergyGraphResult


@pytest.mark.parametrize(
    ("small_molecules", "serialized_molecules", "excluded_residues"),
    [
        (None, None, ["LIG"]),
        ({}, {}, ["LIG"]),
        (
            {"NAD": "NC(=O)c1ccncc1", "ZN": None},
            {"NAD": "NC(=O)c1ccncc1", "ZN": None},
            ["LIG", "NAD"],
        ),
        ({0: "CCO", "HOH": None}, {"0": "CCO", "HOH": None}, ["LIG", 0]),
    ],
)
def test_rbfe_submits_cofactors_and_validates_remaining_residues(
    monkeypatch: pytest.MonkeyPatch,
    small_molecules: dict[str | int, str | None] | None,
    serialized_molecules: dict[str, str | None] | None,
    excluded_residues: list[str | int],
) -> None:
    """Submit fixed cofactors and validate residues without separate parameterization."""
    client = MagicMock()
    client.post.return_value.json.return_value = {
        "name": "test",
        "uuid": "workflow-uuid",
        "created_at": "2026-10-14T00:00:00Z",
        "object_status": 0,
        "parent_uuid": "folder-uuid",
        "notes": "",
        "starred": False,
        "public": False,
        "object_type": "relative_binding_free_energy_perturbation",
        "object_data": {},
        "email_when_complete": False,
        "credits_charged": 0,
        "object_logfile": "",
    }

    @contextmanager
    def mock_api_client() -> Iterator[MagicMock]:
        yield client

    monkeypatch.setattr(rbfe, "api_client", mock_api_client)
    monkeypatch.setattr("rowan.protein.api_client", mock_api_client)
    molecule = rowan.Molecule.from_xyz("2\nhydrogen\nH 0 0 0\nH 0 0 0.74")
    graph_result = RelativeBindingFreeEnergyGraphResult(
        workflow_data={
            "ligands": {"a": molecule.to_stjames(), "b": molecule.to_stjames()},
            "graph": rowan.RBFEGraph(edges=[rowan.RBFEGraphEdge(ligand_a="a", ligand_b="b")]),
        },
        workflow_type="rbfe_graph",
        workflow_uuid="graph-uuid",
    )

    rowan.submit_relative_binding_free_energy_perturbation_workflow(
        graph_result,
        "protein-uuid",
        small_molecules=small_molecules,
    )

    assert client.post.call_args_list[0] == call(
        "/protein/protein-uuid/validate_forcefield", json=excluded_residues
    )
    assert client.post.call_count == 2
    payload = client.post.call_args.kwargs["json"]["workflow_data"]
    assert payload.get("small_molecules") == serialized_molecules
    assert payload["protein"] == "protein-uuid"
    assert set(payload["ligands"]) == {"a", "b"}
