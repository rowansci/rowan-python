"""Tests for macropKa and NMR workflow features."""

from collections.abc import Iterator
from contextlib import contextmanager
from unittest.mock import MagicMock

import stjames
from pytest import MonkeyPatch

import rowan
from rowan.workflows import macropka as macropka_module
from rowan.workflows.nmr import NMRCoupling, NMRResult


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
        "object_type": "macropka",
        "object_data": {},
        "email_when_complete": False,
        "credits_charged": 0,
        "object_logfile": "",
    }


def _mock_submission(monkeypatch: MonkeyPatch) -> MagicMock:
    """Replace the macropKa API client and return its mock client."""
    client = MagicMock()
    client.post.return_value.json.return_value = _workflow_response()

    @contextmanager
    def mock_api_client() -> Iterator[MagicMock]:
        yield client

    monkeypatch.setattr(macropka_module, "api_client", mock_api_client)
    return client


def _molecule() -> stjames.Molecule:
    """Build a molecule with explicit 3D coordinates."""
    return stjames.Molecule.from_xyz("2\nhydrogen\nH 0 0 0\nH 0 0 0.74")


def test_macropka_starling_ii_payload(monkeypatch: MonkeyPatch) -> None:
    """Expose both macroscopic pKa models while retaining Starling as the default."""
    client = _mock_submission(monkeypatch)

    macropka_module.submit_macropka_workflow("CCO", method="starling_ii")
    assert client.post.call_args.kwargs["json"]["workflow_data"]["method"] == "starling_ii"

    macropka_module.submit_macropka_workflow("CCO")
    assert client.post.call_args.kwargs["json"]["workflow_data"]["method"] == "starling"


def test_nmr_couplings_are_typed_and_public() -> None:
    """Return all J-coupling fields through the public result dataclass."""
    workflow = stjames.NMRSpectroscopyWorkflow(
        initial_molecule=_molecule(),
        predicted_couplings=[
            {
                "nuclei": (1, 6),
                "atom_pairs": [(0, 1), (2, 3)],
                "bond_distance": 3,
                "coupling_hz": 7.2,
                "uncertainty_hz": 0.4,
                "conformer_sd_hz": 0.3,
                "model": "coupling-model",
            }
        ],
    )
    result = NMRResult(
        workflow_data=workflow.model_dump(mode="json"),
        workflow_type="nmr",
        workflow_uuid="workflow-uuid",
    )

    assert rowan.NMRCoupling is NMRCoupling
    coupling = result.predicted_couplings[0]
    assert coupling.nuclei == (1, 6)
    assert coupling.atom_pairs == ((0, 1), (2, 3))
    assert coupling.bond_distance == 3
    assert coupling.coupling == 7.2
    assert coupling.uncertainty == 0.4
    assert coupling.conformer_deviation == 0.3
    assert coupling.model == "coupling-model"
