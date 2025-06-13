from typing import Optional

import stjames

from .utils import api_client


class Folder:
    def __init__(
        self,
        uuid: stjames.UUID,
        name: Optional[str] = None,
        parent_uuid: Optional[stjames.UUID] = None,
        notes: Optional[str] = "",
        starred: Optional[bool] = False,
        public: Optional[bool] = False,
    ):
        """
        Initialize a Folder object.

        :param uuid: Unique identifier for the folder.
        :type uuid: UUID
        :param name: Name of the folder.
        :type name: str or None
        :param parent_uuid: UUID of the parent folder.
        :type parent_uuid: UUID or None
        :param notes: Optional notes about the folder.
        :type notes: str
        :param starred: Whether the folder is starred.
        :type starred: bool
        :param public: Whether the folder is public.
        :type public: bool
        """
        self.uuid = uuid
        self.name = name
        self.parent_uuid = parent_uuid
        self.notes = notes
        self.starred = starred
        self.public = public

    def retrieve(self) -> "Folder":
        """
        Retrieve the folder from the API.

        This method refreshes the folder object with the latest data from the API.

        :return: The updated folder object.
        :rtype: Folder
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
    ) -> "Folder":
        """
        Update a folder.

        :param name: The new name of the folder.
        :type name: str or None
        :param parent_uuid: The UUID of the new parent folder.
        :type parent_uuid: stjames.UUID or None
        :param notes: A description of the folder.
        :type notes: str or None
        :param starred: Whether the folder is starred.
        :type starred: bool or None
        :param public: Whether the folder is public.
        :type public: bool or None

        :return: The updated folder object.
        :rtype: Folder
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
    parent_uuid: Optional[stjames.UUID] = None,
    name_contains: Optional[str] = None,
    public: Optional[bool] = None,
    starred: Optional[bool] = None,
    page: int = 0,
    size: int = 10,
) -> list["Folder"]:
    """
    Retrieve a list of folders based on the specified criteria.

    :param parent_uuid: UUID of the parent folder to filter by.
    :param name_contains: Substring to search for in folder names.
    :param public: Filter folders by their public status.
    :param starred: Filter folders by their starred status.
    :param page: Pagination parameter to specify the page number.
    :param size: Pagination parameter to specify the number of items per page.
    :return: A list of Folder objects that match the search criteria.
    :raises: HTTPError if the request to the API fails.
    """

    params = {
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
        items = response.json()["items"]

    return [Folder(**item) for item in items]

def create_folder(
        name: str,
        parent_uuid: Optional[stjames.UUID] = None,
        notes: str = "",
        starred: bool = False,
        public: bool = False,
    ) -> "Folder":
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
