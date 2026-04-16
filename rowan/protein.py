import time
from datetime import datetime
from pathlib import Path
from typing import Any, Self

from pydantic import BaseModel

from .project import Project
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

        :returns: protein with loaded data
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

        :param name: New name of the protein
        :param data: New data of the protein
        :param public: Whether the protein is public
        :param pocket: New pocket of the protein
        :returns: Updated protein object
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

    def sanitize(self, poll_interval: float = 10.0, timeout: float = 300.0) -> None:
        """
        Sanitizes a protein and waits for the process to complete.

        Protein sanitization runs asynchronously on the server. This method
        submits the request then polls until sanitization succeeds, fails, or
        times out.

        :param poll_interval: Seconds between status checks (default 10).
        :param timeout: Maximum seconds to wait before raising (default 300).
        :raises RuntimeError: if sanitization fails, is stopped, or times out.
        :raises requests.HTTPError: if any API request fails.
        """
        with api_client() as client:
            response = client.post(f"/protein/sanitize/{self.uuid}")
            response.raise_for_status()

        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            time.sleep(poll_interval)
            self.refresh()
            match self.sanitized:
                case 2:  # success
                    return
                case 3:  # failed
                    raise RuntimeError(
                        f"Protein sanitization failed for {self.uuid}. "
                        "Check the protein in the Rowan UI for details."
                    )
                case 4:  # stopped
                    raise RuntimeError(f"Protein sanitization was stopped for {self.uuid}.")
                case _:  # 1 (in progress) or None: keep polling
                    pass

        raise RuntimeError(f"Protein sanitization timed out after {timeout:.0f}s for {self.uuid}.")

    def prepare(
        self,
        find_missing_residues: bool = True,
        add_missing_atoms: bool = True,
        remove_heterogens: bool = True,
        keep_waters: bool = False,
        remove_hydrogens: bool = False,
        remove_invalid_hydrogens: bool = False,
        add_hydrogens: bool = True,
        add_hydrogen_ph: float = 7.0,
        optimize_hydrogens: bool = True,
        poll_interval: float = 10.0,
        timeout: float = 300.0,
    ) -> None:
        """
        Prepare a protein for simulation and waits for the process to complete.

        Runs PDBFixer to fix nonstandard residues, add missing atoms/hydrogens,
        and optionally optimizes hydrogen positions with OpenMM. This is the
        recommended method for preparing proteins before MD or RBFE workflows.

        :param find_missing_residues: Identify and model missing residues.
        :param add_missing_atoms: Add missing heavy atoms to residues.
        :param remove_heterogens: Remove ligands, salts, and other heterogens.
        :param keep_waters: Preserve water molecules when removing heterogens.
        :param remove_hydrogens: Remove all existing hydrogens before adding new ones.
        :param remove_invalid_hydrogens: Remove hydrogens not matching the forcefield template.
        :param add_hydrogens: Add missing hydrogen atoms.
        :param add_hydrogen_ph: pH used to determine protonation states when adding hydrogens.
        :param optimize_hydrogens: Optimize hydrogen positions with OpenMM energy minimization.
        :param poll_interval: Seconds between status checks (default 10).
        :param timeout: Maximum seconds to wait before raising (default 300).
        :raises RuntimeError: If preparation fails, is stopped, or times out.
        :raises requests.HTTPError: If any API request fails.
        """
        params = {
            "find_missing_residues": find_missing_residues,
            "add_missing_atoms": add_missing_atoms,
            "remove_heterogens": remove_heterogens,
            "keep_waters": keep_waters,
            "remove_hydrogens": remove_hydrogens,
            "remove_invalid_hydrogens": remove_invalid_hydrogens,
            "add_hydrogens": add_hydrogens,
            "add_hydrogen_ph": add_hydrogen_ph,
            "optimize_hydrogens": optimize_hydrogens,
        }
        with api_client() as client:
            response = client.post(f"/protein/prepare/{self.uuid}", params=params)
            response.raise_for_status()

        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            time.sleep(poll_interval)
            self.refresh()
            match self.sanitized:
                case 2:
                    return
                case 3:
                    raise RuntimeError(
                        f"Protein preparation failed for {self.uuid}. "
                        "Check the protein in the Rowan UI for details."
                    )
                case 4:
                    raise RuntimeError(f"Protein preparation was stopped for {self.uuid}.")
                case _:
                    pass

        raise RuntimeError(f"Protein preparation timed out after {timeout:.0f}s for {self.uuid}.")

    def validate_protein_forcefield(self, exclude_residue_names: list[str] | None = None) -> None:
        """
        Validate that this protein can be parameterized with the MD forcefield.

        Calls the server-side validation which checks that all residues are
        recognized by OpenMM and that there are no clashing atoms. Call this
        before submitting any MD workflow to catch preparation issues early.

        Ligand residues (``LIG``) are always excluded — they are parameterized
        separately by the MD workflow from the provided SMILES.

        If validation fails, try re-preparing with ``remove_invalid_hydrogens=True``:
        ``protein.prepare(remove_invalid_hydrogens=True)``

        :param exclude_residue_names: Additional residue names to skip during validation.
        :raises requests.HTTPError: if validation fails or the API request fails.
        """
        excluded = list({"LIG"} | {name.upper() for name in (exclude_residue_names or [])})
        with api_client() as client:
            response = client.post(
                f"/protein/{self.uuid}/validate_forcefield",
                json=excluded,
            )
            response.raise_for_status()

    def download_pdb_file(self, path: Path | str | None = None, name: str | None = None) -> None:
        """
        Downloads the PDB file for a protein

        :param path: Directory to save the file to (defaults to current directory)
        :param name: Optional custom name for the file (defaults to protein name)
        :raises requests.HTTPError: if the request to the API fails
        """
        path = Path(path) if path is not None else Path.cwd()

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

    :param uuid: UUID of the protein to retrieve.
    :returns: Protein object representing the retrieved protein.
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

    :param ancestor_uuid: UUID of the ancestor protein to filter by
    :param name_contains: Substring to search for in protein names
    :param page: Page number to retrieve
    :param size: Number of items per page
    :returns: List of Protein objects that match the search criteria
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


def upload_protein(
    name: str, file_path: str | Path, project_uuid: str | Project | None = None
) -> Protein:
    """
    Uploads a protein from a PDB file to the API.

    :param name: Name of the protein to create
    :param file_path: Path to the PDB file to upload
    :returns: Protein object representing the uploaded protein
    :raises requests.HTTPError: if the request to the API fails
    """
    file_path = Path(file_path)
    if isinstance(project_uuid, Project):
        project_uuid = project_uuid.uuid
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


def create_protein_from_pdb_id(
    name: str, code: str, project_uuid: str | Project | None = None
) -> Protein:
    """
    Creates a protein from a PDB ID.

    :param name: Name of the protein to create
    :param code: PDB ID of the protein to create
    :param project_uuid: UUID of the project to create the protein in
    :returns: Protein object representing the created protein
    :raises requests.HTTPError: if the request to the API fails
    """
    if isinstance(project_uuid, Project):
        project_uuid = project_uuid.uuid
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
