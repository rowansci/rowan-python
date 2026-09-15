from __future__ import annotations

from datetime import datetime
from typing import TYPE_CHECKING, Any, Self

from pydantic import BaseModel

from .project import default_project, retrieve_project
from .utils import api_client, get_project_uuid

if TYPE_CHECKING:
    from .workflows.base import Workflow


class Folder(BaseModel):
    """A class representing a folder in the Rowan API.

    Attributes:
        uuid: UUID of the folder
        name: name of the folder
        parent_uuid: UUID of the parent folder
        notes: folder notes
        starred: whether the folder is starred
        public: whether the folder is public
        created_at: date and time the folder was created
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
        """Fetch the latest folder data from the API.

        This method refreshes the folder object with the latest data from the API.

        Args:
            in_place: whether to update the current instance in-place

        Returns:
            updated instance (self)

        Raises:
            httpx.HTTPStatusError: API request fails
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
        """Update a folder.

        Args:
            name: new name of the folder
            parent_uuid: UUID of the new parent folder
            notes: description of the folder
            starred: whether the folder is starred
            public: whether the folder is public

        Returns:
            updated folder object
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
        """Delete the folder and all its contents.

        This is a destructive action, it will delete all the folders and
        workflows that are inside this folder.

        Raises:
            httpx.HTTPStatusError: request to the API fails
        """
        with api_client() as client:
            response = client.delete(f"/folder/{self.uuid}")
            response.raise_for_status()

    def print_folder_tree(self, max_depth: int = 10, show_uuids: bool = False) -> None:
        """Retrieves a folder tree from the API.

        Args:
            max_depth: maximum depth of the folder tree
            show_uuids: whether to show the UUIDs of the folders

        Raises:
            httpx.HTTPStatusError: API request fails
        """
        print_folder_tree(self.uuid, max_depth, show_uuids)

    def children(self, size: int = 100) -> list[Folder]:
        """List all child folders directly inside this folder.

        Args:
            size: maximum number of child folders to return

        Returns:
            list of child Folder objects

        Raises:
            httpx.HTTPStatusError: API request fails
        """
        return list_folders(parent_uuid=self.uuid, size=size)

    def workflows(self, size: int = 100) -> list[Workflow]:
        """List all workflows directly inside this folder.

        Args:
            size: maximum number of workflows to return

        Returns:
            list of Workflow objects

        Raises:
            httpx.HTTPStatusError: API request fails
        """
        from .workflows.base import list_workflows

        return list_workflows(parent_uuid=self.uuid, size=size)

    def contents(self, size: int = 100) -> list[Folder | Workflow]:
        """List everything directly inside this folder, both child folders and workflows.

        Folders come first, followed by workflows. For a single type, use `children`
        or `workflows`.

        Args:
            size: maximum number of items of each type to return

        Returns:
            list of Folder and Workflow objects

        Raises:
            httpx.HTTPStatusError: API request fails
        """
        return [*self.children(size=size), *self.workflows(size=size)]

    def parent(self) -> Folder | None:
        """Retrieve the parent folder, or None if this is a root folder.

        Returns:
            parent Folder, or None if there is no parent

        Raises:
            httpx.HTTPStatusError: API request fails
        """
        if self.parent_uuid is None:
            return None
        return retrieve_folder(self.parent_uuid)

    def __truediv__(self, name: str) -> Folder:
        """Traverse into a child folder by name using the `/` operator.

        Examples:
            ```python
            root = rowan.root_folder()
            subfolder = root / "CDK2" / "docking"
            ```

        Args:
            name: exact name of the child folder to navigate into

        Returns:
            child Folder with the given name

        Raises:
            ValueError: no child with that name exists, or if multiple children share
                the same name (use `children` and select by UUID to disambiguate)
        """
        matches = [f for f in self.children(size=200) if f.name == name]
        if not matches:
            raise ValueError(f"No child folder named {name!r} in {self.name!r} ({self.uuid})")
        if len(matches) > 1:
            uuids = ", ".join(f.uuid for f in matches)
            raise ValueError(
                f"Multiple child folders named {name!r} in {self.name!r} ({self.uuid}). "
                f"Use retrieve_folder() with one of these UUIDs to disambiguate: {uuids}"
            )
        return matches[0]


def retrieve_folder(uuid: str) -> Folder:
    """Retrieves a folder from the API by UUID. Folder UUID can be found in the folder's URL.

    Args:
        uuid: UUID of the folder to retrieve

    Returns:
        folder object representing the retrieved folder

    Raises:
        httpx.HTTPStatusError: API request fails
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
    """Retrieve a list of folders based on the specified criteria.

    If no `parent_uuid` is given and a project is active (via `set_project` or
    `rowan.project_uuid`), lists folders rooted at that project's root folder.

    Args:
        parent_uuid: UUID of the parent folder to filter by
        name_contains: substring to search for in folder names
        public: filter folders by their public status
        starred: filter folders by their starred status
        page: pagination parameter to specify the page number
        size: pagination parameter to specify the number of items per page

    Returns:
        list of Folder objects that match the search criteria

    Raises:
        httpx.HTTPStatusError: request to the API fails
    """
    if parent_uuid is None:
        if project_uuid := get_project_uuid():
            parent_uuid = retrieve_project(project_uuid).root_folder_uuid
        else:
            parent_uuid = default_project().root_folder_uuid

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
    """Create a new folder.

    If no `parent_uuid` is given and a project is active (via `set_project` or
    `rowan.project_uuid`), the folder is created inside that project's root folder.

    Args:
        name: name of the folder
        parent_uuid: UUID of the parent folder
        notes: description of the folder
        starred: whether the folder is starred
        public: whether the folder is public

    Returns:
        newly created folder
    """
    if parent_uuid is None:
        if project_uuid := get_project_uuid():
            parent_uuid = retrieve_project(project_uuid).root_folder_uuid
        else:
            parent_uuid = default_project().root_folder_uuid

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


def root_folder() -> Folder:
    """Get the root folder of the active project.

    The root folder is the top of the folder tree you navigate and store workflows in. Use the
    active project set via `set_project` (or `rowan.project_uuid`), falling back to the
    default project.

    Examples:
        ```python
        root = rowan.root_folder()
        for child in root.children():
            print(child.name)
        batch = root / "CDK2" / "docking"
        ```

    Returns:
        root Folder of the active or default project

    Raises:
        httpx.HTTPStatusError: API request fails
    """
    if project_uuid := get_project_uuid():
        root_uuid = retrieve_project(project_uuid).root_folder_uuid
    else:
        root_uuid = default_project().root_folder_uuid
    assert root_uuid is not None
    return retrieve_folder(root_uuid)


def get_folder(path: str, create: bool = True) -> Folder:
    """Get a folder by name or nested path within the default project.

    This is the easiest way to get a folder to use as a location for calculations.
    By default, any missing folders along the path are created automatically.

    Examples:
        ```python
        folder = rowan.get_folder("CDK2/docking/batch_1")
        workflow = rowan.submit_docking_workflow(..., folder_uuid=folder.uuid)
        ```

    Args:
        path: folder name or `/`-separated path, e.g. `"project/subdir/run1"`
        create: create missing folders; otherwise raise ValueError if any
            segment is not found

    Returns:
        deepest `Folder` in the path

    Raises:
        ValueError: path is empty, or `create=False` and a folder is not found
    """
    segments = [s for s in path.split("/") if s]
    if not segments:
        raise ValueError(f"Invalid folder path: {path!r}")
    if project_uuid := get_project_uuid():
        current_uuid = retrieve_project(project_uuid).root_folder_uuid
    else:
        current_uuid = default_project().root_folder_uuid
    folder = None
    for segment in segments:
        matches = list_folders(parent_uuid=current_uuid, name_contains=segment, size=100)
        folder = next(
            (f for f in matches if f.name == segment and f.parent_uuid == current_uuid),
            None,
        )
        if folder is None:
            if not create:
                raise ValueError(f"Folder {segment!r} not found")
            folder = create_folder(name=segment, parent_uuid=current_uuid)
        current_uuid = folder.uuid
    return folder  # type: ignore[return-value]


def print_folder_tree(uuid: str, max_depth: int = 10, show_uuids: bool = False) -> None:
    """Retrieves a folder tree from the API.

    Args:
        uuid: UUID of the root of the folder tree
        max_depth: maximum depth of the folder tree
        show_uuids: whether to show the UUIDs of the folders

    Raises:
        httpx.HTTPStatusError: API request fails
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
