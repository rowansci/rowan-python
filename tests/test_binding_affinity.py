"""Tests for binding affinity workflow submission."""

from collections.abc import Iterator
from contextlib import contextmanager
from typing import Any
from unittest.mock import MagicMock

import pytest
import stjames
from pydantic import ValidationError
from pytest import MonkeyPatch

import rowan
from rowan.workflows import binding_affinity as binding_affinity_module
from rowan.workflows.binding_affinity import BindingAffinityResult, BindingAffinityScore


def _workflow_response() -> dict[str, Any]:
    """Build a minimal binding_affinity workflow API response."""
    return {
        "name": "workflow-1",
        "uuid": "workflow-uuid-1",
        "created_at": "2026-07-17T00:00:00Z",
        "object_status": 0,
        "parent_uuid": "folder-uuid",
        "notes": "",
        "starred": False,
        "public": False,
        "object_type": "binding_affinity",
        "object_data": {},
        "email_when_complete": False,
        "credits_charged": 0,
        "object_logfile": "",
    }


def _mock_client(monkeypatch: MonkeyPatch) -> MagicMock:
    client_mock = MagicMock()
    client_mock.post.return_value.json.return_value = _workflow_response()

    @contextmanager
    def mock_api_client() -> Iterator[MagicMock]:
        yield client_mock

    monkeypatch.setattr(binding_affinity_module, "api_client", mock_api_client)
    return client_mock


@pytest.mark.parametrize(
    ("settings", "settings_json"),
    [
        (rowan.GninaAffinitySettings(), {"settings_type": "gnina"}),
        (rowan.AEVPLIGAffinitySettings(), {"settings_type": "aev_plig"}),
    ],
)
def test_submit_binding_affinity_workflow_holo_protein_ml_settings(
    monkeypatch: MonkeyPatch,
    settings: rowan.BindingAffinitySettings,
    settings_json: dict[str, Any],
) -> None:
    """Submit a holo-protein binding affinity workflow with GNINA/AEV-PLIG settings."""
    client_mock = _mock_client(monkeypatch)

    binding_affinity_module.submit_binding_affinity_workflow(
        protein="protein-uuid",
        ligand_residue_name="LIG",
        binding_affinity_settings=settings,
        folder_uuid="folder-uuid",
    )

    _, kwargs = client_mock.post.call_args
    workflow_data = kwargs["json"]["workflow_data"]
    assert workflow_data["protein"] == "protein-uuid"
    assert workflow_data["ligand_residue_name"] == "LIG"
    assert workflow_data["binding_affinity_settings"] == settings_json
    assert workflow_data["protein_sequences"] == []
    assert workflow_data["ligand_smiles"] == []


def test_submit_binding_affinity_workflow_nesso_sequence_smiles(monkeypatch: MonkeyPatch) -> None:
    """Submit a NESSO binding affinity workflow using sequence/SMILES input, no PDB."""
    client_mock = _mock_client(monkeypatch)

    binding_affinity_module.submit_binding_affinity_workflow(
        protein_sequences=["ACDEFGHIK"],
        ligand_smiles=["CCO"],
        binding_affinity_settings=rowan.NessoAffinitySettings(),
        folder_uuid="folder-uuid",
    )

    _, kwargs = client_mock.post.call_args
    workflow_data = kwargs["json"]["workflow_data"]
    assert workflow_data.get("protein") is None
    assert workflow_data["protein_sequences"] == ["ACDEFGHIK"]
    assert workflow_data["ligand_smiles"] == ["CCO"]
    assert workflow_data["binding_affinity_settings"] == {"settings_type": "nesso"}


def test_submit_binding_affinity_workflow_protein_sequences_requires_nesso() -> None:
    """Reject `protein_sequences` when settings aren't NESSO."""
    with pytest.raises(ValidationError, match=r"protein_sequences.*only supported by NESSO"):
        binding_affinity_module.submit_binding_affinity_workflow(
            protein_sequences=["ACDEFGHIK"],
            ligand_smiles=["CCO"],
            binding_affinity_settings=rowan.GninaAffinitySettings(),
        )


def test_binding_affinity_results_preserve_failed_inputs() -> None:
    """Represent failed rows as `None` without losing input alignment."""
    workflow = stjames.BindingAffinityWorkflow(
        protein_sequences=["ACDEFGHIK"],
        ligand_smiles=["CCO", "CCN"],
        binding_affinity_settings=rowan.NessoAffinitySettings(),
        binding_affinity_results=[{"binding_affinity": -6.5}, None],
    )

    result = BindingAffinityResult(
        workflow_data=workflow.model_dump(mode="json"),
        workflow_type="binding_affinity",
        workflow_uuid="workflow-uuid",
    )

    assert result.scores == [BindingAffinityScore(binding_affinity=-6.5), None]
