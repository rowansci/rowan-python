"""Tests for Rowan request utilities."""

import asyncio

import httpx
from pytest import MonkeyPatch, raises

import rowan
from rowan.utils import (
    api_client,
    api_credentials,
    get_api_key,
    get_project_uuid,
    read_only_api_requests,
)


def test_api_credentials_override_global_configuration(monkeypatch: MonkeyPatch) -> None:
    """Prefer context-local credentials without changing global configuration."""
    monkeypatch.setattr(rowan, "api_key", "global-key")
    monkeypatch.setattr(rowan, "project_uuid", "global-project")
    monkeypatch.setenv("ROWAN_API_KEY", "environment-key")

    with api_credentials("context-key", project_uuid="context-project"):
        assert get_api_key() == "context-key"
        assert get_project_uuid() == "context-project"

    assert get_api_key() == "global-key"
    assert get_project_uuid() == "global-project"


def test_api_credentials_restore_nested_context() -> None:
    """Restore outer credentials after a nested context exits."""
    with api_credentials("outer-key", project_uuid="outer-project"):
        with api_credentials("inner-key"):
            assert get_api_key() == "inner-key"
            assert get_project_uuid() is None

        assert get_api_key() == "outer-key"
        assert get_project_uuid() == "outer-project"


def test_api_credentials_isolate_async_tasks() -> None:
    """Keep credentials isolated across concurrent asynchronous tasks."""

    async def read_after_yield(api_key: str) -> str:
        with api_credentials(api_key):
            await asyncio.sleep(0)
            return get_api_key()

    async def run_workers() -> list[str]:
        return list(await asyncio.gather(read_after_yield("first"), read_after_yield("second")))

    assert asyncio.run(run_workers()) == ["first", "second"]


def test_api_credentials_reject_empty_key() -> None:
    """Reject empty context-local API keys."""
    with raises(ValueError, match="cannot be empty"):
        with api_credentials(""):
            pass


def test_api_credentials_can_hide_the_active_key_from_sdk_helpers() -> None:
    """Allow transports to authenticate without exposing their credential to generic callers."""
    with api_credentials("context-key", reveal_api_key=False):
        with raises(PermissionError, match="not available"):
            get_api_key()


def test_read_only_api_requests_rejects_mutating_http_methods(monkeypatch: MonkeyPatch) -> None:
    """Keep generic integrations from gaining Rowan write authority."""

    def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(200, json={"ok": True})

    original_client = httpx.Client

    def client(*args, **kwargs):
        return original_client(*args, transport=httpx.MockTransport(handler), **kwargs)

    monkeypatch.setattr("rowan.utils.httpx.Client", client)

    with api_credentials("test-key"), read_only_api_requests(), api_client() as http_client:
        assert http_client.get("/workflow").status_code == 200
        with raises(PermissionError, match="only GET/HEAD"):
            http_client.post("/workflow", json={})
