"""Tests for downloading saved RBFE trajectory artifacts."""

from pathlib import Path
from typing import Literal, TypedDict
from unittest.mock import MagicMock

import pytest

from rowan.workflows import relative_binding_free_energy_perturbation as rbfe


class TrajectoryOptions(TypedDict, total=False):
    """Optional trajectory download arguments."""

    leg: Literal["complex", "solvent"]
    lambda_vals: list[float]


@pytest.mark.parametrize(
    ("kwargs", "filename"),
    [
        ({}, "edge_0_trajectories.tar.gz"),
        ({"leg": "solvent", "lambda_vals": [0.0, 1.0]}, "edge_0_solvent_trajectories.tar.gz"),
    ],
)
def test_download_edge_trajectories(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    kwargs: TrajectoryOptions,
    filename: str,
) -> None:
    """Preserve complex defaults and send solvent selections in the API request."""
    expected = tmp_path / filename
    download = MagicMock(return_value=expected)
    monkeypatch.setattr(rbfe, "download_file", download)
    result = MagicMock(
        spec=rbfe.RelativeBindingFreeEnergyPerturbationResult,
        edges=[object()],
        workflow_uuid="workflow-uuid",
    )

    downloaded = rbfe.RelativeBindingFreeEnergyPerturbationResult.download_edge_trajectories(
        result, 0, path=tmp_path, **kwargs
    )

    assert downloaded == expected
    download.assert_called_once_with(
        expected,
        "POST",
        "/trajectory/workflow-uuid/rbfe_trajectory_dcds",
        params={"edge_index": 0, "leg": kwargs.get("leg", "complex")},
        json=kwargs.get("lambda_vals"),
    )
