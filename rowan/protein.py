from typing import Any, Optional

import stjames

from .utils import api_client


class Protein:
    def __init__(
        self,
        uuid: stjames.UUID,
        name: Optional[str] = None,
        data: Optional[dict] = None,
        public: Optional[bool] = None,
        pocket: Optional[list[list[float]]] = None,
    ):
        self.uuid = uuid
        self.name = name
        self.data = data
        self.public = public
        self.pocket = pocket

    def retrieve(self) -> "Protein":
        with api_client() as client:
            response = client.get(f"/protein/{self.uuid}")
            response.raise_for_status()
            protein_data = response.json()

        self.name = protein_data.get("name")
        self.data = protein_data.get("data")
        self.public = protein_data.get("public")
        self.pocket = protein_data.get("pocket")
        return self

    def update(
        self,
        name: Optional[str] = None,
        data: Optional[dict] = None,
        public: Optional[bool] = None,
        pocket: Optional[list[list[float]]] = None,
    ) -> "Protein":
        # Use current values unless new ones are passed in
        updated_payload = {
            "name": name if name is not None else self.name,
            "data": data if data is not None else self.data,
            "public": public if public is not None else self.public,
            "pocket": pocket if pocket is not None else self.pocket,
        }

        with api_client() as client:
            response = client.post(f"/protein/{self.uuid}", json=updated_payload)
            response.raise_for_status()
            updated_data = response.json()

        # Update attributes
        self.name = updated_data.get("name")
        self.data = updated_data.get("data")
        self.public = updated_data.get("public")
        self.pocket = updated_data.get("pocket")
        return self

    def delete(self) -> None:
        with api_client() as client:
            response = client.delete(f"/protein/{self.uuid}")
            response.raise_for_status()

def list_proteins(
    ancestor_uuid: Optional[stjames.UUID] = None,
    name_contains: Optional[str] = None,
    page: int = 0,
    size: int = 20,
) -> list["Protein"]:
    params: dict[str, Any] = {"page": page, "size": size}
    if ancestor_uuid is not None:
        params["ancestor_uuid"] = ancestor_uuid
    if name_contains is not None:
        params["name_contains"] = name_contains

    with api_client() as client:
        response = client.get("/protein", params=params)
        response.raise_for_status()
        results = response.json()["items"]

    return [Protein(**item) for item in results]
