"""Tests for folder path lookup across paginated results."""

from unittest.mock import MagicMock, call

import pytest

import rowan
from rowan import folder as folder_module


@pytest.mark.parametrize("create", [False, True])
def test_get_folder_finds_exact_match_on_later_page(
    monkeypatch: pytest.MonkeyPatch, create: bool
) -> None:
    """Resolve an existing path without raising or creating a duplicate folder."""
    client = MagicMock()
    client.__enter__.return_value = client
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
    monkeypatch.setattr(folder_module, "api_client", lambda: client)
    monkeypatch.setattr(folder_module, "get_project_uuid", lambda: None)
    monkeypatch.setattr(
        folder_module,
        "default_project",
        lambda: rowan.Project(uuid="project", root_folder_uuid="root"),
    )

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


def test_get_folder_creates_only_after_all_pages_are_searched(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Create missing paths after exhausting the matching child folders."""
    client = MagicMock()
    client.__enter__.return_value = client
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
    monkeypatch.setattr(folder_module, "api_client", lambda: client)
    monkeypatch.setattr(folder_module, "get_project_uuid", lambda: None)
    monkeypatch.setattr(
        folder_module,
        "default_project",
        lambda: rowan.Project(uuid="project", root_folder_uuid="root"),
    )

    assert rowan.get_folder("campaign").uuid == "created"
    assert client.get.call_count == 2
    client.post.assert_called_once()
