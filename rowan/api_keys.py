import uuid
from datetime import datetime
from typing import Literal, Self

from pydantic import BaseModel

from .utils import api_client

APIKeyScope = Literal["read", "read_write", "read_write_delete"]


class APIKey(BaseModel):
    """
    Rowan API key.

    :ivar uuid: UUID of the API key.
    :ivar name: Human-readable name of the API key.
    :ivar created_at: When the key was created.
    :ivar expires_at: When the key expires.
    :ivar is_expired: Whether the key has expired.
    :ivar is_revoked: Whether the key has been revoked.
    :ivar scope: Permission scope ("read", "read_write", or "read_write_delete").
    :ivar can_manage_api_keys: Whether this key can create/list/revoke other API keys.
    :ivar budget: Maximum credits this key may spend. ``None`` means unlimited.
    :ivar credits_used: Credits spent by this key so far.
    :ivar scoped_project_uuid: If set, the key can only access this project.
    :ivar created_by_key_uuid: UUID of the API key used to create this one (if any).
    :ivar revoked_at: When the key was revoked, if applicable.
    :ivar last_used_at: When the key was last used, if known.
    """

    uuid: str
    name: str
    created_at: datetime
    expires_at: datetime
    is_expired: bool
    is_revoked: bool
    scope: str
    can_manage_api_keys: bool
    budget: float | None = None
    credits_used: float = 0
    scoped_project_uuid: str | None = None
    created_by_key_uuid: str | None = None
    revoked_at: datetime | None = None
    last_used_at: datetime | None = None

    def __repr__(self) -> str:
        return (
            f"<APIKey name='{self.name}' scope='{self.scope}' "
            f"scoped_project_uuid={self.scoped_project_uuid!r} uuid='{self.uuid}'>"
        )

    def revoke(self) -> Self:
        """
        Revoke this API key.

        :returns: Updated APIKey object.
        """
        with api_client() as client:
            response = client.post(f"/api_key/{self.uuid}/revoke")
            response.raise_for_status()
            return type(self)(**response.json())


class CreatedAPIKey(BaseModel):
    """
    Result of creating a new API key.

    The plaintext ``key`` is only available at creation time — store it now,
    it cannot be retrieved later.

    :ivar key: Plaintext API key. Save this; it is only returned once.
    :ivar api_key: Metadata for the newly-created key.
    """

    key: str
    api_key: APIKey


def create_api_key(
    name: str = "api-key",
    scope: APIKeyScope = "read_write",
    valid_days: int = 365,
    scoped_project_uuid: str | None = None,
    budget: float | None = None,
) -> CreatedAPIKey:
    """
    Create a new API key.

    The caller must currently authenticate with an unscoped key that has
    ``can_manage_api_keys`` permission.

    :param name: Human-readable name for the key.
    :param scope: Permission scope. One of "read", "read_write", "read_write_delete".
    :param valid_days: Number of days until the key expires.
    :param scoped_project_uuid: If provided, restrict the key to a single project.
    :param budget: Maximum credits the key may spend. If not provided, the key has no
        spending limit.
    :returns: plaintext key together with its metadata; the plaintext key is
        only returned once — store it immediately.
    """
    plaintext_key = f"rowan-sk{uuid.uuid4()}"
    payload = {
        "api_key": plaintext_key,
        "name": name,
        "scope": scope,
        "valid_days": valid_days,
        "scoped_project_uuid": scoped_project_uuid,
        "budget": budget,
    }
    with api_client() as client:
        response = client.post("/api_key", json=payload)
        if response.status_code == 403:
            raise PermissionError(
                "API key creation rejected by the server (403). The key you are "
                "authenticating with must be unscoped and have `can_manage_api_keys` "
                "permission. Create a manager key from the Rowan web UI "
                "(Account → API keys) and retry with that key."
            )
        response.raise_for_status()
        return CreatedAPIKey(key=plaintext_key, api_key=APIKey(**response.json()))


def list_api_keys(active: bool | None = True) -> list[APIKey]:
    """
    List API keys belonging to the current user.

    :param active: If True (default), only return non-revoked, non-expired keys.
        If False, only return revoked or expired keys. If None, return all keys.
    :returns: List of APIKey objects.
    """
    params: dict[str, str] = {}
    if active is not None:
        params["active"] = "true" if active else "false"

    with api_client() as client:
        response = client.get("/api_key", params=params)
        response.raise_for_status()
        return [APIKey(**item) for item in response.json()]


def revoke_api_key(uuid: str) -> APIKey:
    """
    Revoke an API key by UUID.

    :param uuid: UUID of the key to revoke.
    :returns: Updated APIKey object.
    """
    with api_client() as client:
        response = client.post(f"/api_key/{uuid}/revoke")
        response.raise_for_status()
        return APIKey(**response.json())


def update_api_key_budget(uuid: str, budget: float | None) -> APIKey:
    """
    Update the spending budget of an API key.

    The caller must currently authenticate with an unscoped key that has
    ``can_manage_api_keys`` permission.

    :param uuid: UUID of the key to update.
    :param budget: New maximum credits the key may spend. Pass ``None`` to remove
        the spending limit.
    :returns: Updated APIKey object.
    """
    with api_client() as client:
        response = client.patch(f"/api_key/{uuid}", json={"budget": budget})
        if response.status_code == 403:
            raise PermissionError(
                "API key update rejected by the server (403). The key you are "
                "authenticating with must be unscoped and have `can_manage_api_keys` "
                "permission. Create a manager key from the Rowan web UI "
                "(Account → API keys) and retry with that key."
            )
        response.raise_for_status()
        return APIKey(**response.json())
