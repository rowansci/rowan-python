"""Base classes for Rowan workflows."""

import logging
import time
import warnings
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, ClassVar, Self, TypeAlias

import stjames
from pydantic import BaseModel, ConfigDict, Field
from rdkit import Chem

from ..molecule import Molecule as RowanMolecule
from ..utils import api_client

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

RdkitMol: TypeAlias = Chem.rdchem.Mol | Chem.rdchem.RWMol
StJamesMolecule: TypeAlias = stjames.Molecule
MoleculeInput: TypeAlias = dict[str, Any] | RowanMolecule | StJamesMolecule | RdkitMol

# Re-export stjames types for convenience
Mode = stjames.Mode
Solvent = stjames.Solvent
MessageType = stjames.MessageType

SolventInput: TypeAlias = Solvent | str | None


@dataclass(frozen=True, slots=True)
class Message:
    """
    A workflow message (error, warning, or info).

    :param title: Short message title.
    :param body: Full message content.
    :param type: Message type: 'error', 'warning', or 'info'.
    """

    title: str
    body: str
    type: str


def parse_messages(raw_messages: list[stjames.Message] | None) -> list[Message]:
    """Parse stjames Message objects into rowan Message objects."""
    if not raw_messages:
        return []
    return [Message(title=m.title, body=m.body, type=m.type) for m in raw_messages]


class WorkflowError(Exception):
    """Raised when a workflow fails or is stopped."""

    pass


@dataclass(slots=True, repr=False)
class WorkflowResult:
    """
    Base class for workflow results.

    Wraps the raw workflow data dict and parses it into a stjames object
    for typed access to nested data.

    :param workflow_data: Raw data dict from the workflow
    :param workflow_type: Workflow type string
    :param workflow_uuid: UUID of the parent workflow (for API calls)
    """

    workflow_data: dict[str, Any]
    workflow_type: str
    workflow_uuid: str
    _workflow: Any = field(default=None, init=False)
    _cache: dict[str, Any] = field(default_factory=dict, init=False)

    _stjames_class: ClassVar[type | None] = None

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}>"

    def __post_init__(self) -> None:
        """Parse workflow data into stjames object for typed access."""
        if self._stjames_class is not None:
            obj = self._stjames_class.model_validate(self.workflow_data)  # type: ignore[attr-defined]
            object.__setattr__(self, "_workflow", obj)

    @property
    def data(self) -> dict[str, Any]:
        """Raw workflow data dict for fallback access."""
        return self.workflow_data

    def clear_cache(self) -> None:
        """
        Clear all cached data to force re-fetching on next access.

        Use this if you need to refresh lazily-loaded data (e.g., structures,
        calculations) from the API.
        """
        self._cache.clear()


RESULT_REGISTRY: dict[str, type[WorkflowResult]] = {}


def register_result(workflow_type: str) -> Callable[[type[WorkflowResult]], type[WorkflowResult]]:
    """Decorator to register a result class for a workflow type."""

    def decorator(cls: type[WorkflowResult]) -> type[WorkflowResult]:
        RESULT_REGISTRY[workflow_type] = cls
        return cls

    return decorator


def create_result(
    workflow_data: dict[str, Any], workflow_type: str, workflow_uuid: str
) -> WorkflowResult:
    """
    Factory function to create the appropriate result type for a workflow.

    :param workflow_data: Raw data dict from the workflow.
    :param workflow_type: Workflow type string.
    :param workflow_uuid: UUID of the parent workflow.
    :returns: Typed WorkflowResult subclass, or base WorkflowResult if unknown.
    """
    result_class = RESULT_REGISTRY.get(workflow_type, WorkflowResult)
    return result_class(
        workflow_data=workflow_data,
        workflow_type=workflow_type,
        workflow_uuid=workflow_uuid,
    )


class Workflow(BaseModel):
    """A Rowan workflow.

    Don't instantiate this class directly. Instead use one of the submit workflow functions.
    Workflow data is not loaded by default to avoid unnecessary downloads that could impact
    performance. Call `fetch_latest()` to fetch and attach the workflow data.

    :param name: Name of the workflow.
    :param uuid: UUID of the workflow.
    :param created_at: Date and time the workflow was created.
    :param updated_at: Date and time the workflow was last updated.
    :param started_at: Date and time the workflow computation was started.
    :param completed_at: Date and time the workflow was completed.
    :param status: Status of the workflow.
    :param parent_uuid: UUID of the parent folder.
    :param notes: Workflow notes.
    :param starred: Whether the workflow is starred.
    :param public: Whether the workflow is public.
    :param workflow_type: Type of the workflow.
    :param data: Data of the workflow.
    :param email_when_complete: Whether to send an email when the workflow completes.
    :param max_credits: Maximum number of credits to use for the workflow.
    :param elapsed: Elapsed time of the workflow.
    :param credits_charged: Number of credits charged for the workflow.
    :param logfile: Workflow logfile.
    """

    name: str
    uuid: str
    created_at: datetime
    updated_at: datetime | None = None
    started_at: datetime | None = None
    completed_at: datetime | None = None
    status: stjames.Status = Field(alias="object_status")
    parent_uuid: str
    notes: str
    starred: bool
    public: bool
    workflow_type: str = Field(alias="object_type")
    data: dict[str, Any] | None = Field(default=None, alias="object_data")
    email_when_complete: bool
    max_credits: int | None = None
    elapsed: float | None = None
    credits_charged: float
    logfile: str = Field(alias="object_logfile")
    compute_hardware: str | None = None

    model_config = ConfigDict(populate_by_name=True)

    def __repr__(self) -> str:
        status = self.status.name.lower() if self.status else "unknown"
        return f"<Workflow name='{self.name}' status={status} uuid='{self.uuid}'>"

    def __str__(self) -> str:
        status = self.status.name.lower() if self.status else "unknown"
        elapsed = f"{self.elapsed:.2f}s" if self.elapsed is not None else "n/a"
        return f"""\
Workflow:  {self.name}
  UUID:    {self.uuid}
  Type:    {self.workflow_type}
  Status:  {status}
  Elapsed: {elapsed}
  Credits: {self.credits_charged}"""

    def fetch_latest(self, in_place: bool = False) -> Self:
        """
        Loads workflow data from the database and updates the current instance.

        :param in_place: Whether to update the current instance in-place.
        :returns: Updated instance (self).
        :raises HTTPError: If the API request fails.
        """
        with api_client() as client:
            response = client.get(f"/workflow/{self.uuid}")
            response.raise_for_status()
            data = response.json()

        if not in_place:
            return type(self).model_validate(data)

        updated_workflow = self.model_validate(data)

        for field_name in type(self).model_fields:
            setattr(self, field_name, getattr(updated_workflow, field_name))

        self.model_rebuild()

        return self

    def update(
        self,
        name: str | None = None,
        parent_uuid: str | None = None,
        notes: str | None = None,
        starred: bool | None = None,
        email_when_complete: bool | None = None,
        public: bool | None = None,
        in_place: bool = False,
    ) -> Self:
        """
        Updates a workflow in the API with new data.

        :param name: New name for the workflow.
        :param parent_uuid: UUID of the parent folder.
        :param notes: Description of the workflow.
        :param starred: Whether the workflow is starred.
        :param email_when_complete: Whether to send an email when complete.
        :param public: Whether the workflow is public.
        :raises HTTPError: If the API request fails.
        """
        old = self.fetch_latest()
        old_data = {
            "name": old.name,
            "parent_uuid": old.parent_uuid,
            "notes": old.notes,
            "starred": old.starred,
            "email_when_complete": old.email_when_complete,
            "public": old.public,
        }

        new_data = {
            "name": name,
            "parent_uuid": parent_uuid,
            "notes": notes,
            "starred": starred,
            "email_when_complete": email_when_complete,
            "public": public,
        }
        updates = {k: v for k, v in new_data.items() if v is not None}
        merged: dict[str, Any] = old_data | updates

        with api_client() as client:
            response = client.post(f"/workflow/{self.uuid}", json=merged)
            response.raise_for_status()
            data = response.json()

        if not in_place:
            return type(self).model_validate(data)

        updated_workflow = self.model_validate(data)

        for field_name in type(self).model_fields:
            setattr(self, field_name, getattr(updated_workflow, field_name))

            self.model_rebuild()

        return self

    def done(self) -> bool:
        """
        Check if the workflow has finished (success, failure, or stopped).

        Non-blocking check following the concurrent.futures.Future pattern.

        :returns: True if workflow is no longer running.
        """
        status = self.get_status()
        return status in {
            stjames.Status.COMPLETED_OK,
            stjames.Status.FAILED,
            stjames.Status.STOPPED,
        }

    def result(self, poll_interval: int = 5) -> "WorkflowResult":
        """
        Wait for completion and return the typed result.

        Follows the concurrent.futures.Future.result() pattern: blocks until
        the workflow completes, then returns a typed result object.

        :param poll_interval: Seconds between status checks while waiting.
        :returns: WorkflowResult subclass with typed access to results.
        :raises WorkflowError: If the workflow failed or was stopped.
        """
        # Wait for completion
        while not self.done():
            time.sleep(poll_interval)

        # Fetch the latest data
        self.fetch_latest(in_place=True)

        # Check for failure
        if self.status in {stjames.Status.FAILED, stjames.Status.STOPPED}:
            status = self.status.name.lower()
            raise WorkflowError(f"Workflow '{self.name}' {status} (uuid={self.uuid})")

        return create_result(self.data or {}, self.workflow_type, self.uuid)

    def peek(self) -> WorkflowResult:
        """
        Get current results without waiting.

        Returns a typed result from whatever data the server has right now.
        Useful for checking progress on long-running jobs or salvaging
        partial results from failed jobs.

        Never raises WorkflowError - just returns what's available.

        :returns: WorkflowResult with current data (may be partial).
        """
        self.fetch_latest(in_place=True)
        return create_result(self.data or {}, self.workflow_type, self.uuid)

    def wait_for_result(self, poll_interval: int = 5) -> Self:
        """
        Wait for the workflow to finish.

        .. deprecated::
            Use :meth:`result` instead, which waits and returns the typed result.

        :returns: Current instance (self).
        """
        while not self.done():
            time.sleep(poll_interval)
        return self

    def get_status(self) -> stjames.Status:
        """
        Gets the status of the workflow.

        :returns: Status of the workflow, as an instance of stjames.Status.
        """
        return self.fetch_latest().status or stjames.Status.QUEUED

    def is_finished(self) -> bool:
        """
        Check if the workflow is finished.

        .. deprecated::
            Use :meth:`done` instead.

        :returns: True if the workflow status is COMPLETED_OK, FAILED, or STOPPED.
        """
        return self.done()

    def stop(self) -> None:
        """
        Stops a workflow.

        :raises HTTPError: If the API request fails.
        """
        with api_client() as client:
            response = client.post(f"/workflow/{self.uuid}/stop")
            response.raise_for_status()

    def delete(self) -> None:
        """
        Deletes the workflow.

        :raises HTTPError: If the API request fails.
        """
        with api_client() as client:
            response = client.delete(f"/workflow/{self.uuid}")
            response.raise_for_status()

    def delete_data(self) -> None:
        """
        Deletes the workflow data from the API.

        :raises HTTPError: If the API request fails.
        """
        with api_client() as client:
            response = client.delete(f"/workflow/{self.uuid}/delete_workflow_data")
            response.raise_for_status()

    def download_msa_files(
        self, msa_format: stjames.MSAFormat, path: Path | str | None = None
    ) -> None:
        """
        Download MSA files for an MSA workflow.

        .. deprecated::
            Use ``workflow.result().download_files()`` instead.
        """
        warnings.warn(
            "download_msa_files() is deprecated. Use workflow.result().download_files() instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        if self.workflow_type != "msa":
            raise ValueError("This workflow is not an MSA workflow.")

        path = Path(path) if path is not None else Path.cwd()

        path.mkdir(parents=True, exist_ok=True)

        with api_client() as client:
            response = client.get(
                f"/workflow/{self.uuid}/get_msa_files",
                params={"msa_format": msa_format.value},
            )
            response.raise_for_status()

        with open(path / f"{self.name}-msa.tar.gz", "wb") as f:
            f.write(response.content)

    def download_dcd_files(
        self, replicates: list[int], name: str | None = None, path: Path | str | None = None
    ) -> None:
        """
        Downloads DCD trajectory files for specified replicates.

        .. deprecated::
            Use ``workflow.result().download_trajectories()`` instead.

        :param replicates: List of replicate indices to download
        :param name: Optional custom name for the tar.gz file
        :param path: Directory to save the file to
        """
        warnings.warn(
            "download_dcd_files() is deprecated. "
            "Use workflow.result().download_trajectories() instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        if self.workflow_type != "pose_analysis_md":
            raise ValueError("This workflow is not a pose analysis molecular dynamics workflow.")

        path = Path(path) if path is not None else Path.cwd()

        path.mkdir(parents=True, exist_ok=True)

        with api_client() as client:
            response = client.post(f"/trajectory/{self.uuid}/trajectory_dcds", json=replicates)
            response.raise_for_status()

        file_path = path / f"{name or self.name}.tar.gz"
        with open(file_path, "wb") as f:
            f.write(response.content)


def extract_smiles(mol: "str | MoleculeInput") -> str:
    """
    Extract a SMILES string from a molecule input or return the string directly.

    :param mol: SMILES string, or any molecule type (RowanMolecule, stjames.Molecule,
        RDKit Mol, or dict).
    :returns: SMILES string.
    :raises TypeError: If the input type is not supported.
    :raises ValueError: If the molecule has no SMILES associated with it.
    """
    if isinstance(mol, str):
        return mol
    elif isinstance(mol, RowanMolecule):
        smiles = mol.smiles
    elif isinstance(mol, stjames.Molecule):
        smiles = mol.smiles
    elif isinstance(mol, dict):
        smiles = mol.get("smiles")
    elif isinstance(mol, Chem.Mol):
        smiles = Chem.MolToSmiles(mol)
    else:
        raise TypeError(f"Cannot extract SMILES from {type(mol)}")
    if smiles is None:
        raise ValueError(
            "Molecule has no SMILES associated with it. Provide a SMILES string directly."
        )
    return smiles


def molecule_to_stjames(mol: MoleculeInput) -> stjames.Molecule:
    """Convert any molecule input type to a stjames.Molecule."""
    match mol:
        case RowanMolecule():
            return mol._to_stjames()
        case stjames.Molecule():
            return mol
        case Chem.Mol():
            return stjames.Molecule.from_rdkit(mol, cid=0)
        case dict():
            return stjames.Molecule(**mol)
        case _:
            raise TypeError(f"Cannot convert {type(mol)} to stjames.Molecule")


def molecule_to_dict(mol: MoleculeInput) -> dict[str, Any]:
    """
    Convert any molecule input type to a dict for API submission.

    :param mol: Molecule as Molecule, stjames.Molecule, RDKit Mol, or dict.
    :returns: Dict representation suitable for API submission.
    """
    match mol:
        case RowanMolecule():
            return mol._to_stjames().model_dump(mode="json")
        case stjames.Molecule():
            return mol.model_dump(mode="json")
        case Chem.Mol():
            return stjames.Molecule.from_rdkit(mol, cid=0).model_dump(mode="json")
        case dict():
            return mol
        case _:
            raise TypeError(f"Cannot convert {type(mol)} to molecule dict")


def submit_workflow(
    workflow_type: stjames.WORKFLOW_NAME,
    workflow_data: dict[str, Any] | None = None,
    initial_molecule: MoleculeInput | None = None,
    initial_smiles: str | None = None,
    name: str | None = None,
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a workflow to the API.

    :param workflow_type: Type of workflow to submit.
    :param workflow_data: Dictionary containing the data required to run the workflow.
    :param initial_molecule: Molecule object to use as the initial molecule.
    :param initial_smiles: SMILES string to use as the initial molecule.
    :param name: Name for the workflow.
    :param folder_uuid: UUID of the folder to store the workflow in.
    :param max_credits: Maximum number of credits to use for the workflow.
    :returns: Workflow object representing the submitted workflow.
    :raises ValueError: If neither `initial_smiles` nor a valid `initial_molecule` is provided.
    :raises HTTPError: If the API request fails.
    """
    if workflow_type not in stjames.WORKFLOW_MAPPING:
        raise ValueError(
            "Invalid workflow type. Must be one of:\n    " + "\n    ".join(stjames.WORKFLOW_MAPPING)
        )

    data: dict[str, Any] = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": workflow_type,
        "workflow_data": workflow_data,
        "max_credits": max_credits,
    }

    if initial_smiles is not None:
        data["initial_smiles"] = initial_smiles
    elif initial_molecule is not None:
        data["initial_molecule"] = molecule_to_dict(initial_molecule)
    else:
        raise ValueError("You must provide either `initial_smiles` or a valid `initial_molecule`.")

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


def retrieve_workflow(uuid: str) -> Workflow:
    """
    Retrieve a workflow from the API by UUID.

    :param uuid: UUID of the workflow to retrieve.
    :returns: Workflow object with the fetched data.
    :raises requests.HTTPError: If the API request fails.
    """
    with api_client() as client:
        response = client.get(f"/workflow/{uuid}")
        response.raise_for_status()
        data = response.json()

    return Workflow.model_validate(data)


def batch_submit_workflow(
    workflow_type: stjames.WORKFLOW_NAME,
    workflow_data: dict[str, Any] | None = None,
    initial_molecules: list[MoleculeInput] | None = None,
    initial_smileses: list[str] | None = None,
    names: list[str] | None = None,
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> list[Workflow]:
    """
    Submits a batch of workflows to the API.

    Each workflow will be submitted with the same workflow type, workflow data,
    and folder UUID, but with different initial molecules and/or SMILES strings.

    :param workflow_type: Type of workflow to submit.
    :param workflow_data: Dictionary containing the data required to run the workflow.
    :param initial_molecules: Molecule objects to use as initial molecules.
    :param initial_smileses: SMILES strings to use as initial molecules.
    :param names: Names for the submitted workflows.
    :param folder_uuid: UUID of the folder to store the workflows in.
    :param max_credits: Maximum number of credits to use per workflow.
    :returns: List of Workflow objects representing the submitted workflows.
    """
    if names is not None:
        expected = len(initial_smileses or initial_molecules or [])
        if len(names) != expected:
            raise ValueError(
                f"Length of names ({len(names)}) must match number of molecules ({expected})."
            )

    workflows = []

    if initial_smileses is not None:
        for i, smiles in enumerate(initial_smileses):
            name = names[i] if names else None
            workflows.append(
                submit_workflow(
                    workflow_type=workflow_type,
                    workflow_data=workflow_data,
                    initial_smiles=smiles,
                    name=name,
                    folder_uuid=folder_uuid,
                    max_credits=max_credits,
                )
            )
    elif initial_molecules is not None:
        for i, molecule in enumerate(initial_molecules):
            name = names[i] if names else None
            workflows.append(
                submit_workflow(
                    workflow_type=workflow_type,
                    workflow_data=workflow_data,
                    initial_molecule=molecule,
                    name=name,
                    folder_uuid=folder_uuid,
                    max_credits=max_credits,
                )
            )
    else:
        raise ValueError("You must provide either `initial_smileses` or `initial_molecules`.")

    return workflows
