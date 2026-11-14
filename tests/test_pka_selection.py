"""Tests for pKa element selection during workflow submission."""

from typing import Literal
from unittest.mock import MagicMock

import pytest

import rowan
from rowan.workflows import pka as pka_module


@pytest.mark.parametrize(
    ("method", "protonate", "deprotonate", "expected_protonate", "expected_deprotonate"),
    [
        ("gxtb_wagen2026", [], [8], [], [8]),
        ("gxtb_wagen2026", [7], [], [7], []),
        ("gxtb_wagen2026", None, None, [7], [7, 8, 16]),
        ("starling", [], [8], [], [8]),
        ("chemprop_nevolianis2025", None, [8], [], [8]),
    ],
)
def test_pka_preserves_explicit_element_selections(
    monkeypatch: pytest.MonkeyPatch,
    method: Literal["gxtb_wagen2026", "starling", "chemprop_nevolianis2025"],
    protonate: list[int] | None,
    deprotonate: list[int] | None,
    expected_protonate: list[int],
    expected_deprotonate: list[int],
) -> None:
    """Default omitted selections while respecting explicitly empty selections."""
    client = MagicMock()
    client.__enter__.return_value = client
    client.post.return_value.json.return_value = {
        "name": "test",
        "uuid": "workflow",
        "created_at": "2026-11-13T00:00:00Z",
        "object_status": 0,
        "parent_uuid": "root",
        "notes": "",
        "starred": False,
        "public": False,
        "object_type": "pka",
        "email_when_complete": False,
        "credits_charged": 0,
        "object_logfile": "",
    }
    monkeypatch.setattr(pka_module, "api_client", lambda: client)
    molecule = rowan.Molecule.from_smiles("NCCO") if method == "gxtb_wagen2026" else "NCCO"

    rowan.submit_pka_workflow(
        molecule, method=method, protonate_elements=protonate, deprotonate_elements=deprotonate
    )

    data = client.post.call_args.kwargs["json"]["workflow_data"]
    assert data["protonate_elements"] == expected_protonate
    assert data["deprotonate_elements"] == expected_deprotonate
