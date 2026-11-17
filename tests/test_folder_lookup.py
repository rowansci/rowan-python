"""Tests for folder path lookup across paginated results."""

from unittest.mock import MagicMock, call

from pytest import MonkeyPatch, fixture, mark, raises

import rowan
from rowan import folder as folder_module


@fixture
def folder_client(monkeypatch: MonkeyPatch) -> MagicMock:
    """Mock folder requests within the default project."""
    client = MagicMock()
    client.__enter__.return_value = client
    monkeypatch.setattr(folder_module, "api_client", lambda: client)
    monkeypatch.setattr(folder_module, "get_project_uuid", lambda: None)
    monkeypatch.setattr(
        folder_module,
        "default_project",
        lambda: rowan.Project(uuid="project", root_folder_uuid="root"),
    )
    return client


@mark.parametrize("create", [False, True])
def test_get_folder_finds_exact_match_on_later_page(folder_client: MagicMock, create: bool) -> None:
    """Resolve an existing path without raising or creating a duplicate folder."""
    client = folder_client
    existing = {"uuid": "existing", "name": "campaign", "parent_uuid": "root"}
    other_matches = [
        {"uuid": f"other-{i}", "name": f"old-{i}-campaign", "parent_uuid": "root"}
        for i in range(100)
    ]
    client.get.return_value.json.side_effect = [
        {"folders": other_matches},
        {"folders": [existing]},
    ]
    client.post.return_value.json.return_value = existing | {"uuid": "duplicate"}
    found = rowan.get_folder("campaign", create=create)

    assert found.uuid == "existing"
    client.post.assert_not_called()
    assert client.get.call_args_list == [
        call(
            "/folder",
            params={"page": page, "size": 100, "parent_uuid": "root", "name_contains": "campaign"},
        )
        for page in range(2)
    ]


@fixture
def missing_folder_client(folder_client: MagicMock) -> MagicMock:
    """Mock paginated results without an exact folder match."""
    client = folder_client
    client.get.return_value.json.side_effect = [
        {
            "folders": [
                {"uuid": f"other-{i}", "name": f"old-{i}-campaign", "parent_uuid": "root"}
                for i in range(100)
            ]
        },
        {"folders": []},
    ]
    client.post.return_value.json.return_value = {
        "uuid": "created",
        "name": "campaign",
        "parent_uuid": "root",
    }
    return client


def test_get_folder_creates_only_after_all_pages_are_searched(
    missing_folder_client: MagicMock,
) -> None:
    """Search every page before creating a missing path."""
    assert rowan.get_folder("campaign").uuid == "created"
    assert missing_folder_client.get.call_count == 2
    missing_folder_client.post.assert_called_once()


def test_get_folder_raises_only_after_all_pages_are_searched(
    missing_folder_client: MagicMock,
) -> None:
    """Search every page before rejecting a missing path without creating it."""
    with raises(ValueError, match="not found"):
        rowan.get_folder("campaign", create=False)
    assert missing_folder_client.get.call_count == 2
    missing_folder_client.post.assert_not_called()


def test_get_folder_restarts_pagination_for_each_segment(folder_client: MagicMock) -> None:
    """Search nested paths from page zero using each resolved parent."""
    folder_client.get.return_value.json.side_effect = [
        {"folders": [{"uuid": "campaign", "name": "campaign", "parent_uuid": "root"}]},
        {
            "folders": [
                {"uuid": f"other-{i}", "name": "old-run", "parent_uuid": "campaign"}
                for i in range(100)
            ]
        },
        {"folders": [{"uuid": "run", "name": "run", "parent_uuid": "campaign"}]},
    ]

    assert rowan.get_folder("campaign/run", create=False).uuid == "run"
    assert folder_client.get.call_args_list == [
        call(
            "/folder",
            params={"page": page, "size": 100, "parent_uuid": parent, "name_contains": name},
        )
        for page, parent, name in [
            (0, "root", "campaign"),
            (0, "campaign", "run"),
            (1, "campaign", "run"),
        ]
    ]
    folder_client.post.assert_not_called()
