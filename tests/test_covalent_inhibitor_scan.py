"""Covalent surface submission, result parsing, and trajectory retrieval."""

from pathlib import Path
from unittest.mock import MagicMock

from pytest import MonkeyPatch, fixture, mark, raises

import rowan
from rowan.workflows import covalent_inhibitor_scan as scan

CYSTEINE = rowan.ResidueRef("polymer", "A", "A.25")
HISTIDINE = rowan.ResidueRef("polymer", "A", "A.162")
LIGAND = rowan.ResidueRef("non_polymer", "A.300")
HYDROGEN = rowan.ResidueAtomRef(HISTIDINE, "HD1")


@fixture
def settings() -> rowan.CovalentInhibitorScanSettings:
    """Two-coordinate surface settings with a descending bond-distance axis."""
    bond = rowan.ScanSettings(
        coordinate=rowan.DistanceCoordinate(
            atoms=(rowan.ResidueAtomRef(CYSTEINE, "SG"), rowan.ResidueAtomRef(LIGAND, "C27"))
        ),
        start=4.6,
        stop=1.6,
        num=3,
        force_constant=100,
    )
    transfer = rowan.ScanSettings(
        coordinate=rowan.LinearCombinationCoordinate(
            coordinates=[
                rowan.DistanceCoordinate(atoms=(rowan.ResidueAtomRef(HISTIDINE, "ND1"), HYDROGEN)),
                rowan.DistanceCoordinate(atoms=(rowan.ResidueAtomRef(LIGAND, "N35"), HYDROGEN)),
            ],
            coefficients=[1, -1],
        ),
        start=-2,
        stop=2,
        num=3,
    )
    return rowan.CovalentInhibitorScanSettings(
        ligand=LIGAND,
        reactive_region=[CYSTEINE, HISTIDINE, LIGAND],
        scan_coordinates=[bond, transfer],
    )


@fixture
def result(settings: rowan.CovalentInhibitorScanSettings) -> rowan.CovalentInhibitorScanResult:
    """Window output with missing cells and a saved shared topology."""
    return rowan.CovalentInhibitorScanResult(
        {
            "protein": "protein-uuid",
            "reactant_smiles": "CC#N",
            "settings": settings.model_dump(mode="json"),
            "topology_uuid": "topology-uuid",
            "windows": [
                {"index": [0, 2], "centers": [4.6, 2], "trajectory": {"uuid": "first"}},
                {"index": [1, 0], "centers": [3.1, -2]},
                {"index": [2, 1], "centers": [1.6, 0], "trajectory": {"uuid": "last"}},
            ],
            "surface": {
                "bin_edges": [[1, 2, 3], [-1, 0, 1]],
                "free_energies": [0, None, 2, 3],
                "effective_samples": [10, 0, 4, 3],
                "uncertainty": [None, None, 0.2, 0.3],
            },
            "messages": [{"title": "Sampling", "body": "More samples needed", "type": "warning"}],
        },
        "covalent_inhibitor_scan",
        "workflow-uuid",
    )


@mark.parametrize(
    ("protein", "dimensions"),
    [("protein-uuid", 1), (rowan.Protein(uuid="protein-uuid"), 2)],
)
def test_submission_preserves_coordinates_and_metadata(
    monkeypatch: MonkeyPatch,
    settings: rowan.CovalentInhibitorScanSettings,
    protein: str | rowan.Protein,
    dimensions: int,
) -> None:
    """Serialize named atom references and weighted coordinates without dropping metadata."""
    settings = settings.model_copy(
        update={"scan_coordinates": settings.scan_coordinates[:dimensions]}
    )
    client = MagicMock()
    client.post.return_value.json.return_value = {
        "uuid": "workflow-uuid",
        "name": "Surface",
        "created_at": "2026-11-12T00:00:00Z",
        "object_status": 6,
        "object_type": "covalent_inhibitor_scan",
        "object_data": {},
        "parent_uuid": "folder-uuid",
        "notes": "",
        "starred": False,
        "public": False,
        "email_when_complete": False,
        "credits_charged": 0,
        "object_logfile": "",
    }
    api = MagicMock()
    api.return_value.__enter__.return_value = client
    monkeypatch.setattr(scan, "api_client", api)
    workflow = rowan.submit_covalent_inhibitor_scan_workflow(
        protein,
        "CC#N",
        settings,
        name="Surface",
        folder=rowan.Folder(uuid="folder-uuid"),
        max_credits=20,
        webhook_url="https://example.com/completed",
        is_draft=True,
    )
    payload = client.post.call_args.kwargs["json"]
    assert client.post.call_args.args == ("/workflow",)
    assert {key: value for key, value in payload.items() if key != "workflow_data"} == {
        "workflow_type": "covalent_inhibitor_scan",
        "initial_smiles": "CC#N",
        "name": "Surface",
        "folder_uuid": "folder-uuid",
        "max_credits": 20,
        "webhook_url": "https://example.com/completed",
        "is_draft": True,
    }
    assert payload["workflow_data"]["protein"] == "protein-uuid"
    assert payload["workflow_data"]["reactant_smiles"] == "CC#N"
    assert payload["workflow_data"]["settings"] == settings.model_dump(mode="json")
    assert workflow.uuid == "workflow-uuid"


def test_result_preserves_grid_order_and_missing_values(
    result: rowan.CovalentInhibitorScanResult,
) -> None:
    """Parse surface cells separately from window order while preserving absent values."""
    assert rowan.create_result(result.data, "covalent_inhibitor_scan", "workflow-uuid") == result
    assert result.surface is not None
    assert [window.index for window in result.windows] == [[0, 2], [1, 0], [2, 1]]
    assert result.windows[1].trajectory is None
    assert result.surface.free_energies == [0, None, 2, 3]
    assert result.surface.uncertainty == [None, None, 0.2, 0.3]
    assert result.messages[0].body == "More samples needed"


def test_topology_uses_workflow_access_and_cache(
    monkeypatch: MonkeyPatch, result: rowan.CovalentInhibitorScanResult
) -> None:
    """Fetch the saved topology with workflow access and avoid repeated requests."""
    topology = rowan.Protein(uuid="topology-uuid")
    retrieve = MagicMock(return_value=topology)
    monkeypatch.setattr(scan, "retrieve_protein", retrieve)
    assert result.get_topology() == topology
    assert result.get_topology() == topology
    retrieve.assert_called_once_with("topology-uuid", workflow_uuid="workflow-uuid")


def test_trajectory_download_keeps_requested_list_positions(
    monkeypatch: MonkeyPatch, tmp_path: Path, result: rowan.CovalentInhibitorScanResult
) -> None:
    """Send list positions unchanged, including reverse order and repeated selections."""
    files = [tmp_path / "last.dcd", tmp_path / "first.dcd", tmp_path / "last.dcd"]
    download = MagicMock(return_value=files)
    monkeypatch.setattr(scan, "download_presigned_files", download)
    assert result.download_trajectory_files([2, 0, 2], str(tmp_path)) == files
    download.assert_called_once_with(
        tmp_path, "/trajectory/workflow-uuid/trajectory_dcd_urls", json=[2, 0, 2]
    )


@mark.parametrize(
    ("indices", "error"), [([-1], IndexError), ([3], IndexError), ([0, 1], ValueError)]
)
def test_trajectory_download_rejects_unavailable_windows(
    monkeypatch: MonkeyPatch,
    result: rowan.CovalentInhibitorScanResult,
    indices: list[int],
    error: type[Exception],
) -> None:
    """Reject missing trajectories or invalid positions before starting any download."""
    download = MagicMock()
    monkeypatch.setattr(scan, "download_presigned_files", download)
    with raises(error):
        result.download_trajectory_files(indices)
    download.assert_not_called()
