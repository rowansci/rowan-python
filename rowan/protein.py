from datetime import datetime
from pathlib import Path
from typing import Any, Self

from pydantic import BaseModel

from .utils import api_client


class Protein(BaseModel):
    """
    A Rowan protein.

    Data is not loaded by default to avoid unnecessary downloads that could impact performance.
    Call `load_data()` to fetch and attach the protein data to this `Protein` object.

    :ivar uuid: The UUID of the protein
    :ivar created_at: The creation date of the protein
    :ivar used_in_workflow: Whether the protein is used in a workflow
    :ivar ancestor_uuid: The UUID of the ancestor protein
    :ivar sanitized: Whether the protein is sanitized
    :ivar name: The name of the protein
    :ivar data: The data of the protein
    :ivar public: Whether the protein is public
    """

    uuid: str
    created_at: datetime | None = None
    used_in_workflow: bool | None = None
    ancestor_uuid: str | None = None
    sanitized: int | None = None
    name: str | None = None
    data: dict | None = None
    public: bool | None = None
    pocket: list[list[float]] | None = None

    def __repr__(self):
        return f"<Protein name='{self.name}' created_at='{self.created_at}' uuid='{self.uuid}'>"

    def refresh(self, in_place: bool = True) -> Self:
        """
        Loads protein data

        :return: protein with loaded data
        """
        with api_client() as client:
            response = client.get(f"/protein/{self.uuid}")
            response.raise_for_status()
            protein_data = response.json()
            if not in_place:
                return self.__class__.model_validate(protein_data)

        self.name = protein_data.get("name")
        self.data = protein_data.get("data")
        self.public = protein_data.get("public")
        self.pocket = protein_data.get("pocket")
        self.sanitized = protein_data.get("sanitized")
        self.used_in_workflow = protein_data.get("used_in_workflow")
        return self

    def update(
        self,
        name: str | None = None,
        data: dict | None = None,
        public: bool | None = None,
        pocket: list[list[float]] | None = None,
    ) -> Self:
        # Use current values unless new ones are passed in
        """
        Updates protein data

        :param name: The new name of the protein
        :param data: The new data of the protein
        :param public: Whether the protein is public
        :param pocket: The new pocket of the protein
        :return: The updated protein object
        """
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
        """
        Deletes a protein

        :raises requests.HTTPError: if the request to the API fails
        """
        with api_client() as client:
            response = client.delete(f"/protein/{self.uuid}")
            response.raise_for_status()

    def sanitize(self) -> None:
        """
        Sanitizes a protein

        :raises requests.HTTPError: if the request to the API fails
        """
        with api_client() as client:
            response = client.post(f"/protein/sanitize/{self.uuid}")
            response.raise_for_status()

    def download_pdb_file(self, path: Path | None = None, name: str | None = None) -> None:
        """
        Downloads the PDB file for a protein

        :param path: Directory to save the file to (defaults to current directory)
        :param name: Optional custom name for the file (defaults to protein name)
        :raises requests.HTTPError: if the request to the API fails
        """
        if path is None:
            path = Path.cwd()

        path.mkdir(parents=True, exist_ok=True)

        with api_client() as client:
            response = client.get(f"/protein/{self.uuid}/get_pdb_file")
            response.raise_for_status()

        file_path = path / f"{name or self.name}.pdb"
        with open(file_path, "wb") as f:
            f.write(response.content)


def retrieve_protein(uuid: str) -> Protein:
    """
    Retrieves a protein from the API using its UUID.

    :param uuid: The UUID of the protein to retrieve.
    :return: A Protein object representing the retrieved protein.
    :raises requests.HTTPError: if the request to the API fails.
    """

    with api_client() as client:
        response = client.get(f"/protein/{uuid}")
        response.raise_for_status()
        protein_data = response.json()

    return Protein(**protein_data)


def list_proteins(
    ancestor_uuid: str | None = None,
    name_contains: str | None = None,
    page: int = 0,
    size: int = 20,
) -> list[Protein]:
    """
    List proteins

    :param ancestor_uuid: The UUID of the ancestor protein to filter by
    :param name_contains: Substring to search for in protein names
    :param page: The page number to retrieve
    :param size: The number of items per page
    :return: A list of Protein objects that match the search criteria
    :raises requests.HTTPError: if the request to the API fails
    """
    params: dict[str, Any] = {"page": page, "size": size}
    if ancestor_uuid is not None:
        params["ancestor_uuid"] = ancestor_uuid
    if name_contains is not None:
        params["name_contains"] = name_contains

    with api_client() as client:
        response = client.get("/protein", params=params)
        response.raise_for_status()
        results = response.json()["proteins"]

    return [Protein(**item) for item in results]


def upload_protein(name: str, file_path: Path, project_uuid: str | None = None) -> Protein:
    """
    Uploads a protein from a PDB file to the API.

    :param name: The name of the protein to create
    :param file_path: The path to the PDB file to upload
    :return: A Protein object representing the uploaded protein
    :raises requests.HTTPError: if the request to the API fails
    """
    with api_client() as client:
        # Step 1: Read the file and post it to the conversion endpoint.
        conversion_payload = {"name": name, "text": file_path.read_text()}
        conversion_response = client.post("/convert/pdb_file_to_protein", json=conversion_payload)
        conversion_response.raise_for_status()  # Ensure the request was successful

        # Extract the JSON data from the conversion response.
        protein_data = conversion_response.json()

        # Step 2: Use the converted data to create the final protein object.
        creation_payload = {
            "name": name,
            "protein_data": protein_data,
            "project_uuid": project_uuid,
        }
        final_response = client.post("/protein", json=creation_payload)
        final_response.raise_for_status()

        # Deserialize the final JSON response into a Protein object and return it.
        return Protein(**final_response.json())


def create_protein_from_pdb_id(name: str, code: str, project_uuid: str | None = None) -> Protein:
    """
    Creates a protein from a PDB ID.

    :param name: The name of the protein to create
    :param code: The PDB ID of the protein to create
    :param project_uuid: The UUID of the project to create the protein in
    :return: A Protein object representing the created protein
    :raises requests.HTTPError: if the request to the API fails
    """
    with api_client() as client:
        # Step 1: Read the file and post it to the conversion endpoint.
        conversion_response = client.post(f"/convert/pdb_id_to_protein?pdb_id={code}")
        conversion_response.raise_for_status()  # Ensure the request was successful

        # Extract the JSON data from the conversion response.
        protein_data = conversion_response.json()

        # Step 2: Use the converted data to create the final protein object.
        creation_payload = {
            "name": name,
            "protein_data": protein_data,
            "project_uuid": project_uuid,
        }
        final_response = client.post("/protein", json=creation_payload)
        final_response.raise_for_status()

        # Deserialize the final JSON response into a Protein object and return it.
        return Protein(**final_response.json())
