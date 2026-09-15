import uuid
from datetime import datetime
from typing import Literal, Self

from pydantic import BaseModel

from .utils import api_client

APIKeyScope = Literal["read", "read_write", "read_write_delete"]


class APIKey(BaseModel):
    """Rowan API key.

    Attributes:
        uuid: UUID of the API key
        name: human-readable name of the API key
        created_at: when the key was created
        expires_at: when the key expires
        is_expired: whether the key has expired
        is_revoked: whether the key has been revoked
        scope: permission scope ("read", "read_write", or "read_write_delete")
        can_manage_api_keys: whether this key can create/list/revoke other API keys
        budget: maximum credits this key may spend. `None` means unlimited
        credits_used: credits spent by this key so far
        scoped_project_uuid: project to which access is restricted
        created_by_key_uuid: UUID of the API key used to create this one (if any)
        revoked_at: when the key was revoked, if applicable
        last_used_at: when the key was last used, if known
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
        """Revoke this API key.

        Returns:
            updated APIKey object
        """
        with api_client() as client:
            response = client.post(f"/api_key/{self.uuid}/revoke")
            response.raise_for_status()
            return type(self)(**response.json())

    def refresh(self, in_place: bool = True) -> Self:
        """Reload this key's metadata (e.g. `credits_used`) from the server.

        Args:
            in_place: update this instance in place rather than return a new instance

        Returns:
            updated APIKey object

        Raises:
            ValueError: this key is no longer present in the account's key list
        """
        matches = [key for key in list_api_keys(active=None) if key.uuid == self.uuid]
        if not matches:
            raise ValueError(f"API key {self.uuid!r} not found.")
        updated = matches[0]

        if not in_place:
            return type(self)(**updated.model_dump())

        for field_name in type(self).model_fields:
            setattr(self, field_name, getattr(updated, field_name))

        return self


class CreatedAPIKey(BaseModel):
    """Result of creating a new API key.

    The plaintext `key` is only available at creation time – store it now,
    it cannot be retrieved later.

    Attributes:
        key: plaintext API key. Save this; it is only returned once
        api_key: metadata for the newly-created key
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
    """Create a new API key.

    The caller must currently authenticate with an unscoped key that has
    `can_manage_api_keys` permission.

    Args:
        name: human-readable name for the key
        scope: permission scope. One of "read", "read_write", "read_write_delete"
        valid_days: number of days until the key expires
        scoped_project_uuid: project to which access is restricted
        budget: maximum credits the key may spend. If not provided, the key has no
            spending limit

    Returns:
        plaintext key together with its metadata; the plaintext key is
        only returned once – store it immediately
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
    """List API keys belonging to the current user.

    Args:
        active: filter by key validity: True for non-revoked, non-expired keys,
            False for revoked or expired keys, or None for all keys

    Returns:
        matching API keys
    """
    params: dict[str, str] = {}
    if active is not None:
        params["active"] = "true" if active else "false"

    with api_client() as client:
        response = client.get("/api_key", params=params)
        response.raise_for_status()
        return [APIKey(**item) for item in response.json()]


def revoke_api_key(uuid: str) -> APIKey:
    """Revoke an API key by UUID.

    Args:
        uuid: UUID of the key to revoke

    Returns:
        updated APIKey object
    """
    with api_client() as client:
        response = client.post(f"/api_key/{uuid}/revoke")
        response.raise_for_status()
        return APIKey(**response.json())


def update_api_key_budget(uuid: str, budget: float | None) -> APIKey:
    """Update the spending budget of an API key.

    The caller must currently authenticate with an unscoped key that has
    `can_manage_api_keys` permission.

    Args:
        uuid: UUID of the key to update
        budget: new maximum credits the key may spend. Pass `None` to remove
            the spending limit

    Returns:
        updated APIKey object
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
