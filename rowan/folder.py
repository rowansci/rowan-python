import stjames
from typing import Optional

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
    ) -> dict:
        with api_client() as client:
            response = client.post(
                "/folder",
                json={
                    "name": name,
                    "parent_uuid": parent_uuid,
                    "notes": notes,
                    "starred": starred,
                    "public": public,
                },
            )
            response.raise_for_status()
            return response.json()

    @classmethod
    def retrieve(cls, uuid: stjames.UUID) -> dict:
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

        new_data = {}
        new_data["name"] = name if name is not None else old_data["name"]
        new_data["parent_uuid"] = (
            parent_uuid if parent_uuid is not None else old_data["parent_uuid"]
        )
        new_data["notes"] = notes if notes is not None else old_data["notes"]
        new_data["starred"] = starred if starred is not None else old_data["starred"]
        new_data["public"] = public if public is not None else old_data["public"]

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
    ):
        params = {"page": page, "size": size}

        if parent_uuid is not None:
            params["parent_uuid"] = parent_uuid

        if name_contains is not None:
            params["name_contains"] = name_contains

        if public is not None:
            params["public"] = public

        if starred is not None:
            params["starred"] = starred

        with api_client() as client:
            response = client.get("/folder", params=params)
            response.raise_for_status()
            return response.json()
