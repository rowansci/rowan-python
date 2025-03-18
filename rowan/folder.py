from typing import Any, Optional

import stjames

from .utils import api_client


class Folder:
    @classmethod
    def create(
        cls,
        name: str,
        parent_uuid: Optional[stjames.UUID] = None,
        notes: str = "",
        starred: bool = False,
        public: bool = False,
    ) -> dict[str, Any]:
        data = {
            "name": name,
            "parent_uuid": parent_uuid,
            "notes": notes,
            "starred": starred,
            "public": public,
        }
        with api_client() as client:
            response = client.post("/folder", json=data)
            response.raise_for_status()
            return response.json()

    @classmethod
    def retrieve(cls, uuid: stjames.UUID) -> dict[str, Any]:
        with api_client() as client:
            response = client.get(f"/folder/{uuid}")
            response.raise_for_status()
            return response.json()

    @classmethod
    def update(
        cls,
        uuid: stjames.UUID,
        name: Optional[str] = None,
        parent_uuid: Optional[stjames.UUID] = None,
        notes: Optional[str] = None,
        starred: Optional[bool] = None,
        public: Optional[bool] = None,
    ) -> None:
        old_data = cls.retrieve(uuid)

        new_data = {
            "name": name if name is not None else old_data["name"],
            "parent_uuid": parent_uuid if parent_uuid is not None else old_data["parent_uuid"],
            "notes": notes if notes is not None else old_data["notes"],
            "starred": starred if starred is not None else old_data["starred"],
            "public": public if public is not None else old_data["public"],
        }

        with api_client() as client:
            response = client.post(f"/folder/{uuid}", json=new_data)
            response.raise_for_status()
            return response.json()

    @classmethod
    def delete(cls, uuid: stjames.UUID) -> None:
        with api_client() as client:
            response = client.delete(f"/folder/{uuid}")
            response.raise_for_status()

    @classmethod
    def list(
        cls,
        parent_uuid: Optional[stjames.UUID] = None,
        name_contains: Optional[str] = None,
        public: Optional[bool] = None,
        starred: Optional[bool] = None,
        page: int = 0,
        size: int = 10,
    ) -> dict[str, Any]:
        params = {
            "page": page,
            "size": size,
            "parent_uuid": parent_uuid,
            "name_contains": name_contains,
            "public": public,
            "starred": starred,
        }

        with api_client() as client:
            response = client.get("/folder", params=params)
            response.raise_for_status()
            return response.json()
