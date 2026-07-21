from datetime import datetime
from typing import Any, Self

from pydantic import BaseModel

import rowan

from .utils import api_client


class Project(BaseModel):
    """
    A class representing a project in the Rowan API.

    :ivar uuid: UUID of the project.
    :ivar name: Name of the project.
    :ivar created_at: Date and time the project was created.
    """

    uuid: str
    name: str | None = None
    root_folder_uuid: str | None = None
    created_at: datetime | None = None

    def __repr__(self) -> str:
        return f"<Project name='{self.name}' created_at='{self.created_at}' uuid='{self.uuid}'>"

    def update(
        self,
        name: str | None = None,
    ) -> Self:
        """
        Update a project.

        :param name: New name of the project.

        :returns: Updated project object.
        """
        payload = {
            "name": name if name is not None else self.name,
        }

        with api_client() as client:
            response = client.post(f"/project/{self.uuid}", json=payload)
            response.raise_for_status()
            updated_data = response.json()

        self.name = updated_data.get("name")
        return self

    def delete(self) -> None:
        """
        Delete the project.

        This is a destructive action, it will delete all the folders and
        workflows that are inside this project.

        :raises requests.HTTPError: if the request to the API fails.
        """
        with api_client() as client:
            response = client.delete(f"/project/{self.uuid}")
            response.raise_for_status()


def retrieve_project(uuid: str) -> Project:
    """
    Retrieves a project from the API by UUID. Project UUID can be found in the project's URL.

    :param uuid: UUID of the project to retrieve.
    :returns: Project object representing the retrieved project.
    :raises HTTPError: If the API request fails.
    """
    with api_client() as client:
        response = client.get(f"/project/{uuid}")
        response.raise_for_status()
        return Project(**response.json())


def list_projects(
    name_contains: str | None = None,
    page: int = 0,
    size: int = 10,
) -> list[Project]:
    """
    Retrieve a list of projects based on the specified criteria.

    :param name_contains: Substring to search for in project names.
    :param page: Pagination parameter to specify the page number.
    :param size: Pagination parameter to specify the number of items per page.
    :returns: List of Folder objects that match the search criteria.
    :raises requests.HTTPError: if the request to the API fails.
    """

    params: dict[str, Any] = {
        "page": page,
        "size": size,
    }

    if name_contains is not None:
        params["name_contains"] = name_contains

    with api_client() as client:
        response = client.get("/project", params=params)
        response.raise_for_status()
        items = response.json()["projects"]

    return [Project(**item) for item in items]


def create_project(
    name: str,
) -> Project:
    """
    Create a new project.

    :param name: Name of the project.
    :returns: Newly created project.
    """
    data = {
        "name": name,
    }
    with api_client() as client:
        response = client.post("/project", json=data)
        response.raise_for_status()
        project_data = response.json()
    return Project(**project_data)


def get_project(name: str, create: bool = False) -> Project:
    """
    Get a project by exact name, optionally creating it if it does not exist.

    The project analogue of :func:`get_folder`. Unlike ``get_folder``, ``create``
    defaults to ``False``: a project is a top-level container, so a typo should
    raise rather than silently spawn a new one. Does not change the active project -
    assign ``rowan.project_uuid`` or use :func:`set_project` for that.

    :param name: Exact name of the project.
    :param create: If True, create the project when no exact match exists.
    :returns: Matched (or newly created) project.
    :raises ValueError: If no match is found and ``create`` is False.
    """
    matches = list_projects(name_contains=name, size=100)
    project = next((p for p in matches if p.name == name), None)
    if project is not None:
        return project
    if create:
        return create_project(name)
    raise ValueError(f"Project {name!r} not found")


def set_project(name: str) -> Project:
    """
    Set the active project by name for all subsequent API calls.

    This is equivalent to setting ``rowan.project_uuid`` directly, but lets
    you use a human-readable name instead of a UUID.

    Example::

        rowan.set_project("CDK2 campaign")
        folder = rowan.get_folder("docking/batch_1")

    :param name: Exact name of the project to activate.
    :returns: Matched project.
    :raises ValueError: If no project with that name is found.
    """
    project = get_project(name)
    rowan.project_uuid = project.uuid
    return project


def default_project() -> Project:
    """
    Retrieves the default project from the API.

    :returns: Project object representing the default project.
    :raises HTTPError: If the API request fails.
    """
    with api_client() as client:
        response = client.get("/user/me/default_project")
        response.raise_for_status()
        return Project(**response.json())
