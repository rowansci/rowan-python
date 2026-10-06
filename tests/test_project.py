"""Tests for project retrieval."""

from collections.abc import Iterator
from contextlib import contextmanager
from unittest.mock import MagicMock

import pytest

from rowan import project as project_module


def _mock_client(monkeypatch: pytest.MonkeyPatch, payload: dict | None) -> None:
    client = MagicMock()
    client.get.return_value.json.return_value = payload

    @contextmanager
    def mock_api_client() -> Iterator[MagicMock]:
        yield client

    monkeypatch.setattr(project_module, "api_client", mock_api_client)


def test_default_project_returns_the_project(monkeypatch: pytest.MonkeyPatch) -> None:
    """Return the default project the API reports."""
    _mock_client(monkeypatch, {"uuid": "project-uuid", "name": "Default"})

    assert project_module.default_project().uuid == "project-uuid"


def test_default_project_without_one_raises_a_clear_error(monkeypatch: pytest.MonkeyPatch) -> None:
    """Explain a missing default project instead of failing to build one from null."""
    _mock_client(monkeypatch, None)

    with pytest.raises(ValueError, match="no default project"):
        project_module.default_project()
