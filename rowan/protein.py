from typing import Any, Optional

import stjames

from .utils import api_client


class Protein:
    @classmethod
    def retrieve(cls, uuid: stjames.UUID) -> dict:
        with api_client() as client:
            response = client.get(f"/protein/{uuid}")
            response.raise_for_status()
            return response.json()

    @classmethod
    def update(
        cls,
        uuid: stjames.UUID,
        name: Optional[str] = None,
        data: Optional[dict] = None,
        public: Optional[bool] = None,
        pocket: Optional[list[list[float]]] = None,
    ) -> None:
        old_data = cls.retrieve(uuid)

        new_data = {}
        new_data["name"] = name if name is not None else old_data["name"]
        new_data["data"] = data if data is not None else old_data["data"]
        new_data["public"] = public if public is not None else old_data["public"]
        new_data["pocket"] = pocket if pocket is not None else old_data["pocket"]

        with api_client() as client:
            response = client.post(f"/protein/{uuid}", json=new_data)
            response.raise_for_status()
            return response.json()

    @classmethod
    def delete(cls, uuid: stjames.UUID) -> None:
        with api_client() as client:
            response = client.delete(f"/protein/{uuid}")
            response.raise_for_status()

    @classmethod
    def list(
        cls,
        ancestor_uuid: Optional[stjames.UUID] = None,
        name_contains: Optional[str] = None,
        page: int = 0,
        size: int = 20,
    ):
        params: dict[str, Any] = {"page": page, "size": size}

        if ancestor_uuid is not None:
            params["ancestor_uuid"] = ancestor_uuid

        if name_contains is not None:
            params["name_contains"] = name_contains

        with api_client() as client:
            response = client.get("/protein", params=params)
            response.raise_for_status()
            return response.json()
