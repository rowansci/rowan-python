"""Regression tests for opt-in Mango force-field support."""

import inspect
from collections.abc import Iterator
from contextlib import contextmanager
from unittest.mock import MagicMock

import stjames
from pytest import MonkeyPatch
from stjames.workflows.relative_binding_free_energy_perturbation import (
    RBFEGraph,
    RBFEGraphEdge,
)

from rowan.workflows import pose_analysis_md as pose_analysis_md_module
from rowan.workflows import protein_md as protein_md_module
from rowan.workflows import relative_binding_free_energy_perturbation as rbfe_module


def _workflow_response(workflow_type: str) -> dict[str, object]:
    """Build a minimal workflow API response."""
    return {
        "name": "test",
        "uuid": "workflow-uuid",
        "created_at": "2026-08-27T00:00:00Z",
        "object_status": 0,
        "parent_uuid": "folder-uuid",
        "notes": "",
        "starred": False,
        "public": False,
        "object_type": workflow_type,
        "object_data": {},
        "email_when_complete": False,
        "credits_charged": 0,
        "object_logfile": "",
    }


def _mock_submission(monkeypatch: MonkeyPatch, module: object, workflow_type: str) -> MagicMock:
    """Replace a workflow module's API client and return its mock client."""
    client = MagicMock()
    client.post.return_value.json.return_value = _workflow_response(workflow_type)

    @contextmanager
    def mock_api_client() -> Iterator[MagicMock]:
        yield client

    monkeypatch.setattr(module, "api_client", mock_api_client)
    return client


def test_protein_md_accepts_mango_without_changing_default(monkeypatch: MonkeyPatch) -> None:
    """Serialize Mango for protein MD while retaining Sage 2.3.0 by default."""
    client = _mock_submission(monkeypatch, protein_md_module, "protein_md")
    protein_md_module.submit_protein_md_workflow(
        "protein-uuid",
        small_molecule_ff="mango_1_0_0",
        validate_forcefield=False,
    )
    payload = client.post.call_args.kwargs["json"]["workflow_data"]
    assert payload["small_molecule_ff"] == "mango_1_0_0"
    assert (
        inspect.signature(protein_md_module.submit_protein_md_workflow)
        .parameters["small_molecule_ff"]
        .default
        == "off_sage_2_3_0"
    )


def test_pose_analysis_md_accepts_mango_without_changing_default(
    monkeypatch: MonkeyPatch,
) -> None:
    """Serialize Mango for pose-analysis MD while retaining Sage 2.3.0 by default."""
    client = _mock_submission(monkeypatch, pose_analysis_md_module, "pose_analysis_md")
    pose_analysis_md_module.submit_pose_analysis_md_workflow(
        "protein-uuid",
        "CC",
        small_molecule_ff="mango_1_0_0",
        validate_forcefield=False,
    )
    payload = client.post.call_args.kwargs["json"]["workflow_data"]
    assert payload["small_molecule_ff"] == "mango_1_0_0"
    assert (
        inspect.signature(pose_analysis_md_module.submit_pose_analysis_md_workflow)
        .parameters["small_molecule_ff"]
        .default
        == "off_sage_2_3_0"
    )


def test_rbfe_accepts_mango_without_changing_default(monkeypatch: MonkeyPatch) -> None:
    """Serialize Mango for RBFE while retaining Sage 2.0.0 by default."""
    client = _mock_submission(
        monkeypatch,
        rbfe_module,
        "relative_binding_free_energy_perturbation",
    )
    molecule = stjames.Molecule.from_xyz("2\nhydrogen\nH 0 0 0\nH 0 0 0.74")
    graph_result = MagicMock()
    graph_result.graph = RBFEGraph(
        edges=[RBFEGraphEdge(ligand_a="ligand_a", ligand_b="ligand_b", core=[(0, 0)])]
    )
    graph_result.ligands = {"ligand_a": molecule, "ligand_b": molecule}

    rbfe_module.submit_relative_binding_free_energy_perturbation_workflow(
        graph_result,
        "protein-uuid",
        tmd_settings="fast",
        forcefield="mango_1_0_0",
        validate_forcefield=False,
    )
    payload = client.post.call_args.kwargs["json"]["workflow_data"]
    assert payload["settings"]["forcefield"] == "mango_1_0_0"
    assert payload["settings"]["charge_method"] == "nagl"
    assert (
        inspect.signature(rbfe_module.submit_relative_binding_free_energy_perturbation_workflow)
        .parameters["forcefield"]
        .default
        == "off_sage_2_0_0"
    )
