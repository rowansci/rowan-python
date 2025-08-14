from datetime import datetime
from typing import Any, Self

from pydantic import BaseModel

from .utils import api_client


class Project(BaseModel):
    """
    A class representing a project in the Rowan API.

    :ivar uuid: The UUID of the project.
    :ivar name: The name of the project.
    :ivar created_at: The date and time the project was created.
    """

    uuid: str
    name: str | None = None
    created_at: datetime | None = None

    def __repr__(self) -> str:
        return f"<Project name='{self.name}' created_at='{self.created_at}' uuid='{self.uuid}'>"

    def update(
        self,
        name: str | None = None,
    ) -> Self:
        """
        Update a project.

        :param name: The new name of the project.

        :return: The updated project object.
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

    :param uuid: The UUID of the project to retrieve.
    :return: A Project object representing the retrieved project.
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
    :return: A list of Folder objects that match the search criteria.
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
        items = response.json()

    return [Project(**item) for item in items]


def create_project(
    name: str,
) -> Project:
    """
    Create a new project.

    :param name: The name of the project.
    :return: The newly created project.
    """
    data = {
        "name": name,
    }
    with api_client() as client:
        response = client.post("/project", json=data)
        response.raise_for_status()
        project_data = response.json()
    return Project(**project_data)


def default_project() -> Project:
    """
    Retrieves the default project from the API.

    :return: A Project object representing the default project.
    :raises HTTPError: If the API request fails.
    """
    with api_client() as client:
        response = client.get("/user/me/default_project")
        response.raise_for_status()
        return Project(**response.json())
