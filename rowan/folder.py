from datetime import datetime
from typing import Any, Self

from pydantic import BaseModel

from .utils import api_client


class Folder(BaseModel):
    """
    A class representing a folder in the Rowan API.

    :ivar uuid: The UUID of the folder.
    :ivar name: The name of the folder.
    :ivar parent_uuid: The UUID of the parent folder.
    :ivar notes: Folder notes.
    :ivar starred: Whether the folder is starred.
    :ivar public: Whether the folder is public.
    :ivar created_at: The date and time the folder was created.
    """

    uuid: str
    name: str | None = None
    parent_uuid: str | None = None
    notes: str = ""
    starred: bool = False
    public: bool = False
    created_at: datetime | None = None

    def __repr__(self) -> str:
        return f"<Folder name='{self.name}' created_at='{self.created_at}' uuid='{self.uuid}'>"

    def fetch_latest(self, in_place: bool = False) -> Self:
        """
        Fetch the latest folder data from the API.

        This method refreshes the folder object with the latest data from the API.

        :param in_place: Whether to update the current instance in-place.
        :return: The updated instance (self).
        :raises HTTPError: If the API request fails.
        """
        with api_client() as client:
            response = client.get(f"/folder/{self.uuid}")
            response.raise_for_status()
            data = response.json()

        if not in_place:
            return self.__class__.model_validate(data)

        updated_folder = self.model_validate(data)

        # Update current instance with new data using class-level model_fields
        for field_name in self.__class__.model_fields:
            setattr(self, field_name, getattr(updated_folder, field_name))

        self.model_rebuild()

        return self

    def update(
        self,
        name: str | None = None,
        parent_uuid: str | None = None,
        notes: str | None = None,
        starred: bool | None = None,
        public: bool | None = None,
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

        :raises requests.HTTPError: if the request to the API fails.
        """
        with api_client() as client:
            response = client.delete(f"/folder/{self.uuid}")
            response.raise_for_status()

    def print_folder_tree(self, max_depth: int = 10, show_uuids: bool = False) -> None:
        """
        Retrieves a folder tree from the API.

        :param max_depth: The maximum depth of the folder tree.
        :param show_uuids: Whether to show the UUIDs of the folders.
        :raises HTTPError: If the API request fails.
        """
        print_folder_tree(self.uuid, max_depth, show_uuids)


def retrieve_folder(uuid: str) -> Folder:
    """
    Retrieves a folder from the API by UUID. Folder UUID can be found in the folder's URL.

    :param uuid: The UUID of the folder to retrieve.
    :return: A Folder object representing the retrieved folder.
    :raises HTTPError: If the API request fails.
    """
    with api_client() as client:
        response = client.get(f"/folder/{uuid}")
        response.raise_for_status()
        return Folder(**response.json())


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
    parent_uuid: str | None = None,
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


def print_folder_tree(uuid: str, max_depth: int = 10, show_uuids: bool = False) -> None:
    """
    Retrieves a folder tree from the API.

    :param uuid: The UUID of the root of the folder tree.
    :param max_depth: The maximum depth of the folder tree.
    :param show_uuids: Whether to show the UUIDs of the folders.
    :raises HTTPError: If the API request fails.
    """
    params: dict[str, Any] = {
        "max_depth": max_depth,
        "show_uuids": show_uuids,
    }
    with api_client() as client:
        response = client.get(f"/folder/{uuid}/folder_tree", params=params)
        response.raise_for_status()
        folder_data = response.json()
    print(folder_data)
