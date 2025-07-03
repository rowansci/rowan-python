from datetime import datetime
from typing import Any, Optional, Self

import stjames
from pydantic import BaseModel

from .utils import api_client


class Folder(BaseModel):
    uuid: stjames.UUID
    name: str | None = None
    parent_uuid: stjames.UUID | None = None
    notes: str = ""
    starred: bool = False
    public: bool = False
    created_at: datetime | None = None

    def __repr__(self) -> str:
        return f"<Folder name='{self.name}' created_at='{self.created_at}'>"

    def retrieve(self) -> Self:
        """
        Retrieve the folder from the API.

        This method refreshes the folder object with the latest data from the API.

        :return: The updated folder object.
        """
        with api_client() as client:
            response = client.get(f"/folder/{self.uuid}")
            response.raise_for_status()
            data = response.json()

        self.name = data.get("name")
        self.parent_uuid = data.get("parent_uuid")
        self.notes = data.get("notes")
        self.starred = data.get("starred")
        self.public = data.get("public")
        return self

    def update(
        self,
        name: Optional[str] = None,
        parent_uuid: Optional[stjames.UUID] = None,
        notes: Optional[str] = None,
        starred: Optional[bool] = None,
        public: Optional[bool] = None,
    ) -> Self:
        """
        Update a folder.

        :param name: The new name of the folder.
        :param parent_uuid: The UUID of the new parent folder.
        :param notes: A description of the folder.
        :param starred: Whether the folder is starred.
        :param public: Whether the folder is public.

        :return: The updated folder object.
        """
        payload = {
            "name": name if name is not None else self.name,
            "parent_uuid": parent_uuid if parent_uuid is not None else self.parent_uuid,
            "notes": notes if notes is not None else self.notes,
            "starred": starred if starred is not None else self.starred,
            "public": public if public is not None else self.public,
        }

        with api_client() as client:
            response = client.post(f"/folder/{self.uuid}", json=payload)
            response.raise_for_status()
            updated_data = response.json()

        self.name = updated_data.get("name")
        self.parent_uuid = updated_data.get("parent_uuid")
        self.notes = updated_data.get("notes")
        self.starred = updated_data.get("starred")
        self.public = updated_data.get("public")
        return self

    def delete(self) -> None:
        """
        Delete the folder and all its contents.

        This is a destructive action, it will delete all the folders and
        workflows that are inside this folder.

        Raises
        ------
            requests.HTTPError: If the request was not successful.
        """
        with api_client() as client:
            response = client.delete(f"/folder/{self.uuid}")
            response.raise_for_status()


def list_folders(
    parent_uuid: str | None = None,
    name_contains: str | None = None,
    public: bool | None = None,
    starred: bool | None = None,
    page: int = 0,
    size: int = 10,
) -> list[Folder]:
    """
    Retrieve a list of folders based on the specified criteria.

    :param parent_uuid: UUID of the parent folder to filter by.
    :param name_contains: Substring to search for in folder names.
    :param public: Filter folders by their public status.
    :param starred: Filter folders by their starred status.
    :param page: Pagination parameter to specify the page number.
    :param size: Pagination parameter to specify the number of items per page.
    :return: A list of Folder objects that match the search criteria.
    :raises requests.HTTPError: if the request to the API fails.
    """

    params: dict[str, Any] = {
        "page": page,
        "size": size,
    }

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
        items = response.json()["folders"]

    return [Folder(**item) for item in items]


def create_folder(
    name: str,
    parent_uuid: Optional[stjames.UUID] = None,
    notes: str = "",
    starred: bool = False,
    public: bool = False,
) -> Folder:
    """
    Create a new folder.

    :param name: The name of the folder.
    :param parent_uuid: The UUID of the parent folder.
    :param notes: A description of the folder.
    :param starred: Whether the folder is starred.
    :param public: Whether the folder is public.
    :return: The newly created folder.
    """
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
        folder_data = response.json()
    return Folder(**folder_data)
