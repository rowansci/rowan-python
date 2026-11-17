"""Tests for basic-calculation workflow submission."""

from collections.abc import Generator
from contextlib import contextmanager
from typing import Any
from unittest.mock import MagicMock

import pytest

import rowan
import stjames
from rowan.workflows import basic_calculation as basic_calculation_module
from rowan.workflows.basic_calculation import PresetName

H2 = stjames.Molecule.from_xyz("H 0 0 0\nH 0 0 0.74")


def _mock_submission(monkeypatch: pytest.MonkeyPatch) -> MagicMock:
    """Replace the basic-calculation API client and return its mock client."""
    client = MagicMock()
    client.post.return_value.json.return_value = {
        "name": "test",
        "uuid": "workflow-uuid",
        "created_at": "2026-11-11T00:00:00Z",
        "object_status": 0,
        "parent_uuid": "folder-uuid",
        "notes": "",
        "starred": False,
        "public": False,
        "object_type": "basic_calculation",
        "object_data": {},
        "email_when_complete": False,
        "credits_charged": 0,
        "object_logfile": "",
    }

    @contextmanager
    def mock_api_client() -> Generator[MagicMock]:
        yield client

    monkeypatch.setattr(basic_calculation_module, "api_client", mock_api_client)
    return client


@pytest.mark.parametrize(
    "thermochem_settings",
    [
        {"temperature": 310, "cutoff_frequency": 50},
        rowan.ThermochemistrySettings(temperature=310, cutoff_frequency=50),
    ],
)
@pytest.mark.parametrize(
    ("method", "engine", "preset"),
    [("GFN2-xTB", "xtb", None), (None, None, "general_nnp")],
)
def test_thermochem_settings_reach_the_workflow(
    monkeypatch: pytest.MonkeyPatch,
    thermochem_settings: dict[str, Any] | rowan.ThermochemistrySettings,
    method: str | None,
    engine: str | None,
    preset: PresetName | None,
) -> None:
    """Thermochemistry settings are submitted with and without a preset."""
    client = _mock_submission(monkeypatch)
    rowan.submit_basic_calculation_workflow(
        H2,
        tasks=["frequencies"],
        thermochem_settings=thermochem_settings,
        method=method,
        engine=engine,
        preset=preset,
    )

    workflow_data = client.post.call_args.kwargs["json"]["workflow_data"]
    assert workflow_data["settings"]["thermochem_settings"] == {
        "temperature": 310,
        "cutoff_frequency": 50,
        "scaling_factor": 1,
    }
