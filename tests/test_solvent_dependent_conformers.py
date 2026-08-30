"""Tests for solvent-dependent conformer workflow updates."""

from collections.abc import Iterator
from contextlib import contextmanager
from unittest.mock import MagicMock

import stjames
from pytest import MonkeyPatch

import rowan
from rowan.workflows import solvent_dependent_conformers as sdc_module
from rowan.workflows.solvent_dependent_conformers import SolventDependentConformersResult


def _workflow_response() -> dict[str, object]:
    """Build a minimal workflow API response."""
    return {
        "name": "test",
        "uuid": "workflow-uuid",
        "created_at": "2026-08-29T00:00:00Z",
        "object_status": 0,
        "parent_uuid": "folder-uuid",
        "notes": "",
        "starred": False,
        "public": False,
        "object_type": "solvent_dependent_conformers",
        "object_data": {},
        "email_when_complete": False,
        "credits_charged": 0,
        "object_logfile": "",
    }


def _mock_submission(monkeypatch: MonkeyPatch) -> MagicMock:
    """Replace the workflow API client and return its mock client."""
    client = MagicMock()
    client.post.return_value.json.return_value = _workflow_response()

    @contextmanager
    def mock_api_client() -> Iterator[MagicMock]:
        yield client

    monkeypatch.setattr(sdc_module, "api_client", mock_api_client)
    return client


def test_submit_solvent_dependent_conformers_tautomer_enumeration(
    monkeypatch: MonkeyPatch,
) -> None:
    """Serialize the default and enabled tautomer-enumeration settings."""
    client = _mock_submission(monkeypatch)
    molecule = rowan.Molecule.from_xyz("2\nhydrogen\nH 0 0 0\nH 0 0 0.74")

    sdc_module.submit_solvent_dependent_conformers_workflow(molecule)
    workflow_data = client.post.call_args.kwargs["json"]["workflow_data"]
    assert workflow_data["enumerate_tautomers"] is False

    sdc_module.submit_solvent_dependent_conformers_workflow(
        molecule,
        enumerate_tautomers=True,
    )
    workflow_data = client.post.call_args.kwargs["json"]["workflow_data"]
    assert workflow_data["enumerate_tautomers"] is True


def test_solvent_dependent_conformer_exposes_tautomer_smiles() -> None:
    """Preserve tautomer identity in conformer results."""
    assert rowan.SMILES is str
    workflow = stjames.SolventDependentConformersWorkflow(
        initial_molecule=stjames.Molecule.from_xyz("2\nhydrogen\nH 0 0 0\nH 0 0 0.74"),
        conformers=[
            {
                "calculation": "calculation-uuid",
                "smiles": "CC(=O)C",
                "free_energy_by_solvent": {"water": -10.0},
                "relative_free_energy_by_solvent": {"water": 0.0},
                "population_by_solvent": {"water": 1.0},
            }
        ],
    )
    result = SolventDependentConformersResult(
        workflow_data=workflow.model_dump(mode="json", serialize_as_any=True),
        workflow_type="solvent_dependent_conformers",
        workflow_uuid="workflow-uuid",
    )

    assert result.conformers[0].smiles == "CC(=O)C"
