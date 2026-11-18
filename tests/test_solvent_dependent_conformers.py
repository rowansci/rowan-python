"""Tests for solvent-dependent conformer workflow updates."""

from collections.abc import Generator
from contextlib import contextmanager
from unittest.mock import MagicMock

import pytest

import rowan
import stjames
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


def _mock_submission(monkeypatch: pytest.MonkeyPatch) -> MagicMock:
    """Replace the workflow API client and return its mock client."""
    client = MagicMock()
    client.post.return_value.json.return_value = _workflow_response()

    @contextmanager
    def mock_api_client() -> Generator[MagicMock]:
        yield client

    monkeypatch.setattr(sdc_module, "api_client", mock_api_client)
    return client


def test_submit_solvent_dependent_conformers_tautomer_enumeration(
    monkeypatch: pytest.MonkeyPatch,
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
    assert result.conformers[0].singlepoint_energy is None
    assert result.conformers[0].rrho_correction is None
    assert result.conformers[0].solvation_free_energy_by_solvent == {}


def test_conformer_energy_precision_and_components() -> None:
    """Expose microhartree precision and the decomposition used for solvent scoring."""
    workflow = stjames.SolventDependentConformersWorkflow(
        initial_molecule=stjames.Molecule.from_xyz("2\nhydrogen\nH 0 0 0\nH 0 0 0.74"),
        conformers=[
            {
                "calculation": "calculation-uuid",
                "free_energy_by_solvent": {"water": -100.11790015},
                "relative_free_energy_by_solvent": {"water": 0},
                "population_by_solvent": {"water": 1},
                "singlepoint_energy": -100.1234567,
                "rrho_correction": 0.01234567,
                "solvation_free_energy_by_solvent": {"water": -0.00678912},
            }
        ],
    )
    result = SolventDependentConformersResult(
        workflow.model_dump(mode="json"), "solvent_dependent_conformers", "workflow-uuid"
    )
    conformer = result.conformers[0]
    assert conformer.free_energy_by_solvent == {rowan.Solvent.WATER: -100.117900}
    assert conformer.singlepoint_energy == -100.123457
    assert conformer.rrho_correction == 0.012346
    assert conformer.solvation_free_energy_by_solvent == {rowan.Solvent.WATER: -0.006789}
