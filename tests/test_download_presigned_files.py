"""Tests for downloading files through presigned links."""

from contextlib import contextmanager
from pathlib import Path
from typing import Any, Generator

import httpx
import pytest

from rowan import utils


def test_download_presigned_files_refreshes_expired_links(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Fetch links from the API, stream them without the API key, and refresh on 403."""
    link_batches = iter(
        [
            [
                {"filename": "a.dcd", "url": "https://s3/a?1"},
                {"filename": "../b.dcd", "url": "https://s3/b?1"},
            ],
            [
                {"filename": "a.dcd", "url": "https://s3/a?2"},
                {"filename": "../b.dcd", "url": "https://s3/b?2"},
            ],
        ]
    )
    api_posts: list[tuple[str, Any]] = []

    class FakeAPIClient:
        def post(self, url: str, params: Any = None, json: Any = None) -> httpx.Response:
            api_posts.append((url, json))
            return httpx.Response(200, json=next(link_batches))

    @contextmanager
    def fake_api_client() -> Generator[FakeAPIClient, None, None]:
        yield FakeAPIClient()

    s3_requests: list[httpx.Request] = []

    def s3(request: httpx.Request) -> httpx.Response:
        s3_requests.append(request)
        if str(request.url) == "https://s3/b?1":
            return httpx.Response(403)
        return httpx.Response(200, content=request.url.path.encode())

    real_client = httpx.Client
    monkeypatch.setattr(utils, "api_client", fake_api_client)
    monkeypatch.setattr(
        utils.httpx,
        "Client",
        lambda **kwargs: real_client(transport=httpx.MockTransport(s3), **kwargs),
    )

    paths = utils.download_presigned_files(
        tmp_path / "out", "/trajectory/x/trajectory_dcd_urls", json=[0, 1]
    )

    assert paths == [tmp_path / "out" / "a.dcd", tmp_path / "out" / "b.dcd"]
    assert [p.read_bytes() for p in paths] == [b"/a", b"/b"]
    assert len(api_posts) == 2
    assert [str(r.url) for r in s3_requests] == [
        "https://s3/a?1",
        "https://s3/b?1",
        "https://s3/b?2",
    ]
    assert all("x-api-key" not in r.headers for r in s3_requests)
