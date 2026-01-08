import logging
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Literal, Self, TypeAlias

import stjames
from pydantic import BaseModel, Field
from rdkit import Chem
from stjames.optimization.freezing_string_method import FSMSettings

from .protein import Protein
from .utils import api_client

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

RdkitMol: TypeAlias = Chem.rdchem.Mol | Chem.rdchem.RWMol
StJamesMolecule: TypeAlias = stjames.Molecule


class Workflow(BaseModel):
    """A Rowan workflow.

    Don’t instantiate this class directly. Instead use one of the submit workflow functions.
    Workflow data is not loaded by default to avoid unnecessary downloads that could impact
    performance. Call `load_data()` to fetch and attach the workflow data to this `Workflow` object.

    :ivar name: The name of the workflow.
    :var uuid: The UUID of the workflow.
    :var created_at: The date and time the workflow was created.
    :var updated_at: The date and time the workflow was last updated.
    :var started_at: The date and time the workflow computation was started.
    :var completed_at: The date and time the workflow was completed.
    :var status: The status of the workflow.
    :var parent_uuid: The UUID of the parent folder.
    :var notes: Workflow notes.
    :var starred: Whether the workflow is starred.
    :var public: Whether the workflow is public.
    :var workflow_type: The type of the workflow.
    :var data: The data of the workflow.
    :var email_when_complete: Whether the workflow should send an email when it is complete.
    :var max_credits: The maximum number of credits to use for the workflow.
    :var elapsed: The elapsed time of the workflow.
    :var credits_charged: The number of credits charged for the workflow.
    :var logfile: The workflow's logfile.
    ...
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

    class Config:  # noqa: D106
        validate_by_name = True

    def __repr__(self) -> str:
        return f"<Workflow name='{self.name}' created_at='{self.created_at}' uuid='{self.uuid}'>"

    def fetch_latest(self, in_place: bool = False) -> Self:
        """
        Loads workflow data from the database and updates the current instance.

        :param in_place: Whether to update the current instance in-place.
        :return: The updated instance (self).
        :raises HTTPError: If the API request fails.
        """
        with api_client() as client:
            response = client.get(f"/workflow/{self.uuid}")
            response.raise_for_status()
            data = response.json()

        if not in_place:
            return self.__class__.model_validate(data)

        # Create a new instance with proper field mapping
        updated_workflow = self.model_validate(data)

        # Update current instance with new data using class-level model_fields
        for field_name in self.__class__.model_fields:
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

        :param name: The new name of the workflow.
        :param parent_uuid: The UUID of the parent folder.
        :param notes: A description of the workflow.
        :param starred: Whether the workflow is starred.
        :param email_when_complete: Whether the workflow should send an email when it is complete.
        :param public: Whether the workflow is public.
        :raises HTTPError: If the API request fails.
        """
        old_data = self.fetch_latest()

        new_data = {
            "name": name if name is not None else old_data.name,
            "parent_uuid": parent_uuid if parent_uuid is not None else old_data.parent_uuid,
            "notes": notes if notes is not None else old_data.notes,
            "starred": starred if starred is not None else old_data.starred,
            "email_when_complete": email_when_complete
            if email_when_complete is not None
            else old_data.email_when_complete,
            "public": public if public is not None else old_data.public,
        }

        with api_client() as client:
            response = client.post(f"/workflow/{self.uuid}", json=new_data)
            response.raise_for_status()
            data = response.json()

        if not in_place:
            return self.__class__.model_validate(data)

        updated_workflow = self.model_validate(data)

        for field_name in self.__class__.model_fields:
            setattr(self, field_name, getattr(updated_workflow, field_name))

            self.model_rebuild()

        return self

    def wait_for_result(self, poll_interval: int = 5) -> Self:
        """
        Wait for the workflow to finish.

        This method will block until the workflow has finished computing. It
        will periodically poll the API to check the status of the workflow.

        :return: The current instance, with the workflow data loaded.
        """

        if poll_interval <= 0:
            raise ValueError("poll_interval must be a positive integer")

        while not self.is_finished():
            time.sleep(poll_interval)
        return self

    def get_status(self) -> stjames.Status:
        """
        Gets the status of the workflow.

        :return: The status of the workflow, as an instance of stjames.Status.
        """
        return self.fetch_latest().status or stjames.Status.QUEUED

    def is_finished(self) -> bool:
        """
        Check if the workflow is finished.

        This method checks the current status of the workflow and determines if it has completed,
        failed, or been stopped.

        :return: True if the workflow status is COMPLETED_OK, FAILED, or STOPPED; False otherwise.
        """

        status = self.get_status()

        return status in {
            stjames.Status.COMPLETED_OK,
            stjames.Status.FAILED,
            stjames.Status.STOPPED,
        }

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

    def download_msa_files(self, msa_format: stjames.MSAFormat, path: Path | None = None) -> None:
        if self.workflow_type != "msa":
            raise ValueError("This workflow is not an MSA workflow.")

        if path is None:
            path = Path.cwd()

        if not path.exists():
            path.mkdir(parents=True, exist_ok=True)

        with api_client() as client:
            response = client.get(
                f"/workflow/{self.uuid}/get_msa_files", params={"msa_format": msa_format.value}
            )
            response.raise_for_status()

        with open(path / f"{self.name}-msa.tar.gz", "wb") as f:
            f.write(response.content)

    def download_dcd_files(
        self, replicates: list[int], name: str | None = None, path: Path | None = None
    ) -> None:
        """
        Downloads DCD trajectory files for specified replicates

        :param replicates: List of replicate indices to download
        :param name: Optional custom name for the tar.gz file (defaults to workflow name)
        :param path: Directory to save the file to (defaults to current directory)
        :raises ValueError: if workflow is not a pose analysis MD workflow
        :raises requests.HTTPError: if the request to the API fails
        """
        if self.workflow_type != "pose_analysis_md":
            raise ValueError("This workflow is not a pose analysis molecular dynamics workflow.")

        if path is None:
            path = Path.cwd()

        path.mkdir(parents=True, exist_ok=True)

        with api_client() as client:
            response = client.post(f"/trajectory/{self.uuid}/trajectory_dcds", json=replicates)
            response.raise_for_status()

        file_path = path / f"{name or self.name}.tar.gz"
        with open(file_path, "wb") as f:
            f.write(response.content)


def submit_workflow(
    workflow_type: stjames.WORKFLOW_NAME,
    workflow_data: dict[str, Any] | None = None,
    initial_molecule: dict[str, Any] | StJamesMolecule | RdkitMol | None = None,
    initial_smiles: str | None = None,
    name: str | None = None,
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a workflow to the API.

    :param workflow_type: The type of workflow to submit.
    :param workflow_data: A dictionary containing the data required to run the workflow.
    :param initial_molecule: A molecule object to use as the initial molecule in the workflow.
    At least one of a molecule or SMILES must be provided.
    :param initial_smiles: A SMILES string to use as the initial molecule in the workflow.
    At least one of a molecule or SMILES must be provided.
    :param name: A name to give to the workflow.
    :param folder_uuid: The UUID of the folder to store the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises ValueError: If neither `initial_smiles` nor a valid `initial_molecule` is provided.
    :raises HTTPError: If the API request fails.
    """
    if workflow_type not in stjames.WORKFLOW_MAPPING:
        raise ValueError(
            "Invalid workflow type. Must be one of:\n    " + "\n    ".join(stjames.WORKFLOW_MAPPING)
        )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": workflow_type,
        "workflow_data": workflow_data,
        "max_credits": max_credits,
    }

    if initial_smiles is not None:
        data["initial_smiles"] = initial_smiles
    elif isinstance(initial_molecule, StJamesMolecule):
        data["initial_molecule"] = initial_molecule.model_dump(mode="json")
    elif isinstance(initial_molecule, dict):
        data["initial_molecule"] = initial_molecule
    elif isinstance(initial_molecule, RdkitMol):
        data["initial_molecule"] = StJamesMolecule.from_rdkit(initial_molecule, cid=0).model_dump(
            mode="json"
        )
    else:
        raise ValueError("You must provide either `initial_smiles` or a valid `initial_molecule`.")

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


def batch_submit_workflow(
    workflow_type: stjames.WORKFLOW_NAME,
    workflow_data: dict[str, Any] | None = None,
    initial_molecules: list[dict[str, Any] | StJamesMolecule | RdkitMol] | None = None,
    initial_smileses: list[str] | None = None,
    names: list[str] | None = None,
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> list[Workflow]:
    """
    Submits a batch of workflows to the API. Each workflow will be submitted with the
    same workflow type, workflow data, and folder UUID, but with different initial molecules and/or
    SMILES strings.

    :param workflow_type: The type of workflow to submit.
    :param workflow_data: A dictionary containing the data required to run the workflow.
    :param initial_molecules: A list of molecule objects to use as the initial molecules in the
    workflows.
    At least one of a molecule or SMILES must be provided.
    :param initial_smileses: A list of SMILES strings to use as the initial molecules in the
    workflows.
    At least one of a molecule or SMILES must be provided.
    :param names: A list of names to give to the workflows.
    :param folder_uuid: The UUID of the folder to store the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow. This applies per
    workflow, not for the batch.
    :return: A list of Workflow objects representing the submitted workflows.
    :raises ValueError: If neither `initial_smiles` nor a valid `initial_molecule` is provided.
    :raises HTTPError: If the API request fails.
    """
    if workflow_type not in stjames.WORKFLOW_MAPPING:
        raise ValueError(
            "Invalid workflow type. Must be one of:\n    " + "\n    ".join(stjames.WORKFLOW_MAPPING)
        )

    if initial_molecules and initial_smileses:
        raise ValueError(
            "You must provide either `initial_molecules` or `initial_smileses`, but not both."
        )

    if names:
        if initial_molecules and len(names) != len(initial_molecules):
            logger.warning(
                "The number of names is not the same as the number of initial "
                "molecules. Generic names of the form 'Batch Workflow {i}' "
                "will be used."
            )
        if initial_smileses and len(names) != len(initial_smileses):
            logger.warning(
                "The number of names is not the same as the number of initial SMILES."
                "Generic names of the form 'Batch Workflow {i}' will be used."
            )

    data = {
        "names": names,
        "folder_uuid": folder_uuid,
        "workflow_type": workflow_type,
        "workflow_data": workflow_data,
        "max_credits": max_credits,
        "initial_molecules": initial_molecules,
        "initial_smileses": initial_smileses,
    }

    with api_client() as client:
        response = client.post("/workflow/batch_submit", json=data)
        response.raise_for_status()
        return [Workflow(**workflow_data) for workflow_data in response.json()]


def retrieve_workflow(uuid: str) -> Workflow:
    """
    Retrieves a workflow from the API.

    :param uuid: The UUID of the workflow to retrieve.
    :return: A Workflow object representing the retrieved workflow.
    :raises HTTPError: If the API request fails.
    """
    with api_client() as client:
        response = client.get(f"/workflow/{uuid}")
        response.raise_for_status()
        return Workflow(**response.json())


def retrieve_workflows(uuids: list[str]) -> list[Workflow]:
    """
    Retrieves a list of workflows from the API.

    :param uuids: The UUIDs of the workflows to retrieve.
    :return: A list of Workflow objects representing the retrieved workflows.
    :raises HTTPError: If the API request fails.
    """
    with api_client() as client:
        response = client.post("/workflow/batch_retrieve", json={"uuids": uuids})
        response.raise_for_status()
        return [Workflow(**workflow_data) for workflow_data in response.json()]


def batch_poll_status(uuids: list[str]) -> list[Workflow]:
    """
    Polls the status of a list of workflows.

    :param uuids: The UUIDs of the workflows to poll.
    :return: A dictionary of statuses and the count of workflows in that status.
    :raises HTTPError: If the API request fails.
    """
    with api_client() as client:
        response = client.post("/workflow/batch_status", json={"uuids": uuids})
        response.raise_for_status()
        return response.json()


def retrieve_calculation_molecules(
    uuid: str, return_frequencies: bool = False
) -> list[dict[str, Any]]:
    """
    Retrieves a list of molecules from a calculation.

    :param uuid: The UUID of the calculation to retrieve molecules from.
    :param return_frequencies: Whether to return the frequencies of the molecules.
    :return: A list of dictionaries representing the molecules in the calculation.
    :raises HTTPError: If the API request fails.
    """
    with api_client() as client:
        response = client.get(
            f"/calculation/{uuid}/molecules", params={"return_frequencies": return_frequencies}
        )
        response.raise_for_status()
        return response.json()


def list_workflows(
    parent_uuid: str | None = None,
    name_contains: str | None = None,
    public: bool | None = None,
    starred: bool | None = None,
    status: int | None = None,
    workflow_type: stjames.WORKFLOW_NAME | None = None,
    page: int = 0,
    size: int = 10,
) -> list[Workflow]:
    """
    Lists workflows subject to the specified criteria.

    :param parent_uuid: The UUID of the parent folder.
    :param name_contains: Substring to search for in workflow names.
    :param public: Filter workflows by their public status.
    :param starred: Filter workflows by their starred status.
    :param status: Filter workflows by their status.
    :param workflow_type: Filter workflows by their type.
    :param page: The page number to retrieve.
    :param size: The number of items per page.
    :return: A list of Workflow objects that match the search criteria.
    :raises requests.HTTPError: if the request to the API fails.
    """
    params: dict[str, Any] = {"page": page, "size": size}

    if parent_uuid is not None:
        params["parent_uuid"] = parent_uuid

    if name_contains is not None:
        params["name_contains"] = name_contains

    if public is not None:
        params["public"] = public

    if starred is not None:
        params["starred"] = starred

    if status is not None:
        params["object_status"] = status

    if workflow_type is not None:
        params["object_type"] = workflow_type

    with api_client() as client:
        response = client.get("/workflow", params=params)
        response.raise_for_status()
        return [Workflow(**item) for item in response.json()["workflows"]]


def submit_basic_calculation_workflow(
    initial_molecule: dict[str, Any] | StJamesMolecule | RdkitMol,
    method: stjames.Method | str = "uma_m_omol",
    basis_set: stjames.BasisSet | str | None = None,
    tasks: list[str] | None = None,
    mode: str = "auto",
    engine: str | None = None,
    name: str = "Basic Calculation Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submit a basic calculation workflow to the API.

    :param initial_molecule: The molecule to perform the calculation on.
    :param method: The method to use for the calculation.
    :param basis_set: The basis_set to use (if any).
    for options.
    :param tasks: A list of tasks to perform for the calculation.
    :param mode: The mode to run the calculation in. See [list of available modes](https://github.com/rowansci/stjames-public/blob/master/stjames/mode.py)
    for options.
    :param engine: The engine to use for the calculation. See [list of available engines](https://github.com/rowansci/stjames-public/blob/master/stjames/engine.py)
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if not tasks:
        tasks = ["optimize"]

    if isinstance(initial_molecule, StJamesMolecule):
        initial_molecule = initial_molecule.model_dump(mode="json")
    elif isinstance(initial_molecule, RdkitMol):
        initial_molecule = StJamesMolecule.from_rdkit(initial_molecule, cid=0)

    if isinstance(method, str):
        method = stjames.Method(method)

    workflow = stjames.BasicCalculationWorkflow(
        initial_molecule=initial_molecule,
        settings=stjames.Settings(
            method=method,
            basis_set=basis_set,
            tasks=tasks,
            mode=mode,
        ),
        engine=engine or method.default_engine(),
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "basic_calculation",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": initial_molecule,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


def submit_conformer_search_workflow(
    initial_molecule: dict[str, Any] | StJamesMolecule | RdkitMol,
    conf_gen_mode: str = "rapid",
    final_method: stjames.Method | str = "aimnet2_wb97md3",
    solvent: str | None = None,
    transition_state: bool = False,
    name: str = "Conformer Search Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a conformer search workflow to the API.

    :param initial_molecule: The molecule to perform the conformer search on.
    :param conf_gen_mode: The mode to use for conformer generation. See
    [list of available modes](https://github.com/rowansci/stjames-public/blob/master/stjames/mode.py)
    for options.
    :param final_method: The method to use for the final optimization.
    See [list of available methods](https://github.com/rowansci/stjames-public/blob/master/stjames/method.py)
    for options.
    :param solvent: The solvent to use for the final optimization. See [the list of available solvents](https://github.com/rowansci/stjames-public/blob/master/stjames/solvent.py)
        for valid values. Be aware that not all methods support solvents.
    :param transition_state: Whether to optimize the transition state.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if isinstance(initial_molecule, StJamesMolecule):
        initial_molecule = initial_molecule.model_dump(mode="json")
    elif isinstance(initial_molecule, RdkitMol):
        initial_molecule = StJamesMolecule.from_rdkit(initial_molecule, cid=0)

    if isinstance(final_method, str):
        final_method = stjames.Method(final_method)

    solvent_model = None
    if solvent:
        solvent_model = "alpb" if final_method in stjames.XTB_METHODS else "cpcm"

    opt_settings = stjames.Settings(
        method=final_method,
        tasks=["optimize"],
        mode=stjames.Mode.AUTO,
        solvent_settings={"solvent": solvent, "model": solvent_model} if solvent else None,
        opt_settings={"transition_state": transition_state, "constraints": []},
    )

    msos = stjames.MultiStageOptSettings(
        mode=stjames.Mode.MANUAL,
        xtb_preopt=True,
        optimization_settings=[opt_settings],
    )

    workflow = stjames.ConformerSearchWorkflow(
        initial_molecule=initial_molecule,
        multistage_opt_settings=msos,
        conf_gen_mode=conf_gen_mode,
        mso_mode="manual",
        solvent=solvent,
        transition_state=transition_state,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "conformer_search",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": initial_molecule,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


def submit_solubility_workflow(
    initial_smiles: str,
    solubility_method: Literal["fastsolv", "kingfisher", "esol"] = "fastsolv",
    solvents: list[str] | None = None,
    temperatures: list[float] | None = None,
    name: str = "Solubility Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a solubility workflow to the API.

    :param solubility_method: The name of the desired model for solubility prediction.
    :param initial_smiles: The smiles of the molecule to calculate the solubility of.
    :param solvents: The list of solvents to use for the calculation.
    :param temperatures: The list of temperatures to use for the calculation.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """

    if not solvents:
        solvents = ["CCCCCC", "CC1=CC=CC=C1", "C1CCCO1", "CC(=O)OCC", "CCO", "CC#N"]

    if not temperatures:
        temperatures = [273.15, 298.15, 323.15, 348.15, 373.15]

    workflow = stjames.SolubilityWorkflow(
        initial_smiles=initial_smiles,
        solubility_method=solubility_method,
        solvents=solvents,
        temperatures=temperatures,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "solubility",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_smiles": initial_smiles,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


def submit_pka_workflow(
    initial_molecule: dict[str, Any] | StJamesMolecule | RdkitMol | str,
    pka_range: tuple[int, int] = (2, 12),
    method: Literal["aimnet2_wagen2024", "chemprop_nevolianis2025"] = "aimnet2_wagen2024",
    solvent: str = "water",
    deprotonate_elements: list[int] | None = None,
    protonate_elements: list[int] | None = None,
    mode: str = "careful",
    name: str = "pKa Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a pKa workflow to the API.

    :param initial_molecule: The molecule to calculate the pKa of.
        Valid input types include `stjames.Molecule`, RDKit molecule, and SMILES.
    :param pka_range: The range of pKa values to calculate.
    :param method: The algorithm used to compute pKa values.
    :param solvent: The solvent in which pKa values will be computed.
    :param deprotonate_elements: The elements to deprotonate. Given by atomic number.
    :param protonate_elements: The elements to protonate. Given by atomic number.
    :param mode: The mode to run the calculation in. See
    [list of available modes](https://github.com/rowansci/stjames-public/blob/master/stjames/mode.py)
    for options.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    initial_smiles: str = ""
    initial_stjames_mol: StJamesMolecule | None = None

    if isinstance(initial_molecule, dict):
        initial_stjames_mol = StJamesMolecule.model_validate(initial_molecule)
    elif isinstance(initial_molecule, StJamesMolecule):
        initial_stjames_mol = initial_molecule
    elif isinstance(initial_molecule, RdkitMol):
        initial_stjames_mol = StJamesMolecule.from_rdkit(initial_molecule, cid=0)
    elif isinstance(initial_molecule, str):
        initial_smiles = initial_molecule

    protonate_elements = protonate_elements or [7]
    deprotonate_elements = deprotonate_elements or [7, 8, 16]

    workflow = stjames.pKaWorkflow(
        initial_molecule=initial_stjames_mol,
        initial_smiles=initial_smiles,
        pka_range=pka_range,
        deprotonate_elements=deprotonate_elements,
        protonate_elements=protonate_elements,
        mode=mode,
        solvent=solvent,
        microscopic_pka_method=method,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "pka",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": initial_stjames_mol.model_dump(mode="json")
        if initial_stjames_mol
        else None,
        "initial_smiles": initial_smiles,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


def submit_redox_potential_workflow(
    initial_molecule: dict[str, Any] | StJamesMolecule | RdkitMol,
    reduction: bool = False,
    oxidization: bool = True,
    mode: str = "rapid",
    name: str = "Redox Potential Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a redox potential workflow to the API.

    :param initial_molecule: The molecule to calculate the redox potential of.
    :param reduction: Whether to calculate the reduction potential.
    :param oxidization: Whether to calculate the oxidization potential.
    :param mode: The mode to run the calculation in. See
    [list of available modes](https://github.com/rowansci/stjames-public/blob/master/stjames/mode.py)
    for options.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if isinstance(initial_molecule, StJamesMolecule):
        initial_molecule = initial_molecule.model_dump(mode="json")
    elif isinstance(initial_molecule, RdkitMol):
        initial_molecule = StJamesMolecule.from_rdkit(initial_molecule, cid=0)

    workflow = stjames.RedoxPotentialWorkflow(
        initial_molecule=initial_molecule,
        oxidation=oxidization,
        reduction=reduction,
        mode=mode,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "redox_potential",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": initial_molecule,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


def submit_fukui_workflow(
    initial_molecule: dict[str, Any] | StJamesMolecule | RdkitMol,
    optimization_method: str = "gfn2_xtb",
    fukui_method: str = "gfn1_xtb",
    solvent_settings: dict[str, str] | None = None,
    name: str = "Fukui Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a fukui workflow to the API.

    :param initial_molecule: The molecule to calculate the fukui indices of.
    :param optimization_method: The method to use for the optimization.
    :param fukui_method: The method to use for the fukui calculation.
    :param solvent_settings: The solvent settings to use for the fukui calculation.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if isinstance(initial_molecule, StJamesMolecule):
        initial_molecule = initial_molecule.model_dump(mode="json")
    elif isinstance(initial_molecule, RdkitMol):
        initial_molecule = StJamesMolecule.from_rdkit(initial_molecule, cid=0)

    optimization_settings = stjames.Settings(method=optimization_method)
    fukui_settings = stjames.Settings(method=fukui_method, solvent_settings=solvent_settings)

    stjames.FukuiIndexWorkflow(
        initial_molecule=initial_molecule,
        optimization_settings=optimization_settings,
        fukui_settings=fukui_settings,
    )

    workflow_data = {
        "opt_settings": optimization_settings.model_dump(mode="json"),
        "opt_engine": stjames.Method(optimization_method).default_engine(),
        "fukui_settings": fukui_settings.model_dump(mode="json"),
        "fukui_engine": stjames.Method(fukui_method).default_engine(),
    }

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "fukui",
        "workflow_data": workflow_data,
        "initial_molecule": initial_molecule,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


def submit_tautomer_search_workflow(
    initial_molecule: dict[str, Any] | StJamesMolecule | RdkitMol,
    mode: str = "careful",
    name: str = "Tautomer Search Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a tautomer search workflow to the API.

    :param initial_molecule: The molecule to calculate the tautomers of.
    :param mode: The mode to run the calculation in. See
    [list of available modes](https://github.com/rowansci/stjames-public/blob/master/stjames/mode.py)
    for options.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if isinstance(initial_molecule, StJamesMolecule):
        initial_molecule = initial_molecule.model_dump(mode="json")
    elif isinstance(initial_molecule, RdkitMol):
        initial_molecule = StJamesMolecule.from_rdkit(initial_molecule, cid=0)

    workflow = stjames.TautomerWorkflow(
        initial_molecule=initial_molecule,
        mode=mode,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "tautomers",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": initial_molecule,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


def submit_descriptors_workflow(
    initial_molecule: dict[str, Any] | StJamesMolecule | RdkitMol,
    name: str = "Descriptors Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a descriptors workflow to the API.

    :param initial_molecule: The molecule to calculate the descriptors of.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if isinstance(initial_molecule, StJamesMolecule):
        initial_molecule = initial_molecule.model_dump(mode="json")
    elif isinstance(initial_molecule, RdkitMol):
        initial_molecule = StJamesMolecule.from_rdkit(initial_molecule, cid=0)

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "descriptors",
        "workflow_data": {},
        "initial_molecule": initial_molecule,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


def submit_scan_workflow(
    initial_molecule: dict[str, Any] | StJamesMolecule | RdkitMol,
    scan_settings: stjames.ScanSettings | dict[str, Any] | None = None,
    calculation_engine: str | None = None,
    calculation_method: stjames.Method | str = "uma_m_omol",
    wavefront_propagation: bool = True,
    name: str = "Scan Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a scan workflow to the API.

    :param initial_molecule: The molecule used in the scan.
    :param scan_settings: The scan settings.
    :param calculation_engine: The engine to use for the calculation. See [list of available engines](https://github.com/rowansci/stjames-public/blob/master/stjames/engine.py)
    :param calculation_method: The method to use for the calculation.
    See [list of available methods](https://github.com/rowansci/stjames-public/blob/master/stjames/method.py)
    for options.
    :param wavefront_propagation: Whether to use wavefront propagation in the scan.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to store the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if isinstance(initial_molecule, StJamesMolecule):
        initial_molecule = initial_molecule.model_dump(mode="json")
    elif isinstance(initial_molecule, RdkitMol):
        initial_molecule = StJamesMolecule.from_rdkit(initial_molecule, cid=0)

    if isinstance(calculation_method, str):
        calculation_method = stjames.Method(calculation_method)

    calc_settings = stjames.Settings(
        method=calculation_method,
        tasks=["optimize"],
        corrections=[],
        mode="auto",
        opt_settings={"constraints": []},
    )

    workflow = stjames.ScanWorkflow(
        initial_molecule=initial_molecule,
        scan_settings=scan_settings,
        calc_settings=calc_settings,
        calc_engine=calculation_engine or calculation_method.default_engine(),
        wavefront_propagation=wavefront_propagation,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "scan",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": initial_molecule,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


def submit_macropka_workflow(
    initial_smiles: str,
    min_pH: int = 0,
    max_pH: int = 14,
    min_charge: int = -2,
    max_charge: int = 2,
    compute_solvation_energy: bool = False,
    compute_aqueous_solubility: bool = False,
    name: str = "Macropka Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a macropka workflow to the API.

    :param initial_smiles: The molecule used in the macropka workflow.
    :param min_pH: The minimum pH to use in the macropka workflow.
    :param max_pH: The maximum pH to use in the macropka workflow.
    :param min_charge: The minimum charge to use in the macropka workflow.
    :param max_charge: The maximum charge to use in the macropka workflow.
    :param compute_aqueous_solubility: Whether to compute the aqueous solubility for each pH.
    :param compute_solvation_energy: Whether to compute the solvation energy.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to store the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """

    workflow = stjames.MacropKaWorkflow(
        initial_smiles=initial_smiles,
        min_pH=min_pH,
        max_pH=max_pH,
        min_charge=min_charge,
        max_charge=max_charge,
        compute_solvation_energy=compute_solvation_energy,
        compute_aqueous_solubility=compute_aqueous_solubility,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "macropka",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_smiles": initial_smiles,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


def submit_irc_workflow(
    initial_molecule: dict[str, Any] | StJamesMolecule | RdkitMol | None = None,
    method: stjames.Method | str = "uma_m_omol",
    preopt: bool = True,
    step_size: float = 0.05,
    max_irc_steps: int = 30,
    name: str = "IRC Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits an Intrinsic Reaction Coordinate (IRC) workflow to the API.

    :param initial_molecule: The initial molecule to perform the IRC calculation on.
    :param method: The computational method to use for the IRC calculation.
    See [list of available methods](https://github.com/rowansci/stjames-public/blob/master/stjames/method.py)
    for options.
    :param preopt: Whether to perform a pre-optimization of the molecule.
    :param step_size: The step size to use for the IRC calculation.
    :param max_irc_steps: The maximum number of IRC steps to perform.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted IRC workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """

    if isinstance(initial_molecule, StJamesMolecule):
        initial_molecule = initial_molecule.model_dump(mode="json")
    elif isinstance(initial_molecule, RdkitMol):
        initial_molecule = StJamesMolecule.from_rdkit(initial_molecule, cid=0)

    if isinstance(method, str):
        method = stjames.Method(method)

    workflow = stjames.IRCWorkflow(
        initial_molecule=initial_molecule,
        settings=stjames.Settings(
            method=method,
            tasks=[],
            corrections=[],
            mode="auto",
        ),
        preopt=preopt,
        step_size=step_size,
        max_irc_steps=max_irc_steps,
        mode="manual",
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "irc",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": initial_molecule,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


def submit_protein_cofolding_workflow(
    initial_protein_sequences: list[str],
    initial_smiles_list: list[str] | None = None,
    ligand_binding_affinity_index: int | None = None,
    use_msa_server: bool = True,
    use_potentials: bool = False,
    compute_strain: bool = False,
    do_pose_refinement: bool = False,
    name: str = "Cofolding Workflow",
    model: str = stjames.CofoldingModel.BOLTZ_2.value,
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a protein cofolding workflow to the API.

    :param initial_protein_sequences: The sequences of the proteins to be cofolded.
    :param initial_smiles_list: A list of SMILES strings for the ligands to be cofolded with.
    :param ligand_binding_affinity_index: The index of the ligand for which to compute the binding affinity.
    :param use_msa_server: Whether to use the MSA server for the computation.
    :param use_potentials: Whether to use potentials for the computation.
    :param do_pose_refinement: whether to optimize non-rotatable bonds in output poses
    :param compute_strain: whether to compute the strain of the pose (if `pose_refinement` is enabled)
    :param name: The name of the workflow.
    :param model: The model to use for the computation.
    :param folder_uuid: The UUID of the folder to store the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """  # noqa: E501

    workflow = stjames.ProteinCofoldingWorkflow(
        use_msa_server=use_msa_server,
        use_potentials=use_potentials,
        model=model,
        ligand_binding_affinity_index=ligand_binding_affinity_index,
        initial_smiles_list=initial_smiles_list,
        initial_protein_sequences=initial_protein_sequences,
        do_pose_refinement=do_pose_refinement,
        compute_strain=compute_strain,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "protein_cofolding",
        "workflow_data": workflow.model_dump(mode="json"),
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


def submit_docking_workflow(
    protein: str | Protein,
    pocket: list[list[float]],
    initial_molecule: dict[str, Any] | StJamesMolecule | RdkitMol | None = None,
    executable: str = "vina",
    scoring_function: str = "vinardo",
    exhaustiveness: float = 8,
    do_csearch: bool = False,
    do_optimization: bool = False,
    do_pose_refinement: bool = False,
    name: str = "Docking Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a Docking workflow to the API.

    :param protein: The protein to dock. Can be input as a uuid or a Protein object.
    :param initial_molecule: The initial molecule to be docked
    :param executable: Which docking implementation to use.
    :param scoring_function: Which docking scoring function to use.
    :param exhaustiveness: Which exhaustiveness to employ.
    :param do_csearch: Whether to perform a conformational search on the ligand.
    :param do_optimization: Whether to perform an optimization on the ligand.
    :param do_pose_refinement: Whether or not to optimize output poses.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted docking workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """

    if isinstance(initial_molecule, StJamesMolecule):
        initial_molecule = initial_molecule.model_dump(mode="json")
    elif isinstance(initial_molecule, RdkitMol):
        initial_molecule = StJamesMolecule.from_rdkit(initial_molecule, cid=0)

    if isinstance(protein, Protein):
        protein = protein.uuid

    docking_settings = {
        "executable": executable,
        "exhaustiveness": exhaustiveness,
        "scoring_function": scoring_function,
    }

    workflow = stjames.DockingWorkflow(
        initial_molecule=initial_molecule,
        target_uuid=protein,
        pocket=pocket,
        do_csearch=do_csearch,
        do_optimization=do_optimization,
        do_pose_refinement=do_pose_refinement,
        docking_settings=docking_settings,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "docking",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": initial_molecule,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


def submit_ion_mobility_workflow(
    initial_molecule: dict[str, Any] | StJamesMolecule | RdkitMol,
    temperature: float = 300,
    protonate: bool = False,
    do_csearch: bool = True,
    do_optimization: bool = True,
    name: str = "Ion-Mobility Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits an ion-mobility workflow to the API.

    :param initial_molecule: The molecule used in the scan.
    :param temperature: The temperature at which to predict CCS values.
    :param protonate: Whether or not to automatically detect protonation site.
        If `True`, every basic site will be protonated and values returned for the most stable.
    :param do_csearch: Whether to perform a conformational search on the molecule.
    :param do_optimization: Whether to perform an optimization on the molecule.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to store the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if isinstance(initial_molecule, StJamesMolecule):
        initial_molecule = initial_molecule.model_dump(mode="json")
    elif isinstance(initial_molecule, RdkitMol):
        initial_molecule = StJamesMolecule.from_rdkit(initial_molecule, cid=0)

    workflow = stjames.IonMobilityWorkflow(
        initial_molecule=initial_molecule,
        temperature=temperature,
        protonate=protonate,
        do_csearch=do_csearch,
        do_optimization=do_optimization,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "ion_mobility",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": initial_molecule,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


def submit_nmr_workflow(
    initial_molecule: dict[str, Any] | StJamesMolecule | RdkitMol,
    solvent: str | None = "chloroform",
    do_csearch: bool = True,
    do_optimization: bool = True,
    name: str = "NMR Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits an NMR-prediction workflow to the API.

    :param initial_molecule: The molecule used in the scan.
    :param solvent: The solvent in which to compute NMR spectra.
    :param do_csearch: Whether to perform a conformational search on the input structure.
    :param do_optimization: Whether to perform an optimization on the input structure.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to store the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if isinstance(initial_molecule, StJamesMolecule):
        initial_molecule = initial_molecule.model_dump(mode="json")
    elif isinstance(initial_molecule, RdkitMol):
        initial_molecule = StJamesMolecule.from_rdkit(initial_molecule, cid=0)

    workflow_data = {"initial_molecule": initial_molecule, "solvent": solvent}

    if not do_csearch:
        workflow_data["conf_gen_settings"] = None

    if not do_optimization:
        workflow_data["multistage_opt_settings"] = None

    workflow = stjames.NMRSpectroscopyWorkflow.model_validate(workflow_data)

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "nmr",
        "workflow_data": workflow.model_dump(serialize_as_any=True, mode="json"),
        "initial_molecule": initial_molecule,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


def submit_strain_workflow(
    initial_molecule: dict[str, Any] | StJamesMolecule | RdkitMol,
    name: str = "Strain Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a strain workflow to the API.

    :param initial_molecule: The molecule used in the scan.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to store the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if isinstance(initial_molecule, StJamesMolecule):
        initial_molecule = initial_molecule.model_dump(mode="json")
    elif isinstance(initial_molecule, RdkitMol):
        initial_molecule = StJamesMolecule.from_rdkit(initial_molecule, cid=0)

    workflow = stjames.StrainWorkflow(initial_molecule=initial_molecule)

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "strain",
        "workflow_data": workflow.model_dump(serialize_as_any=True, mode="json"),
        "initial_molecule": initial_molecule,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


def submit_double_ended_ts_search_workflow(
    reactant: dict[str, Any] | StJamesMolecule,
    product: dict[str, Any] | StJamesMolecule,
    calculation_settings: stjames.Settings | dict[str, Any] | None = None,
    search_settings: FSMSettings | dict[str, Any] | None = None,
    optimize_inputs: bool = False,
    optimize_ts: bool = True,
    name: str = "Double-Ended TS Search Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a double-ended transition state search workflow to the API.

    :param reactant: reactant Molecule.
    :param product: product Molecule.
    :param calculation_settings: Settings to use for calculations.
    :param search_settings: settings to use for the transition state search.
    :param optimize_inputs: Whether to optimize the reactant and product before the search.
    :param optimize_ts: Whether to optimize the found transition state.
    :param name: name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: Workflow object representing the submitted workflow.
    """
    workflow = stjames.DoubleEndedTSSearchWorkflow(
        reactant=reactant,
        product=product,
        calculation_settings=calculation_settings,
        search_settings=search_settings,
        optimize_inputs=optimize_inputs,
        optimize_ts=optimize_ts,
    )
    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "double_ended_ts_search",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": reactant
        if isinstance(reactant, dict)
        else reactant.model_dump(mode="json"),
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


def submit_pose_analysis_md_workflow(
    protein: str | Protein,
    initial_smiles: str,
    num_trajectories: int = 1,
    simulation_time_ns: int = 10,
    ligand_residue_name: str = "LIG",
    name: str = "Pose-Analysis MD Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a pose-analysis MD workflow to the API.

    :param protein: The *holo* protein on which MD will be run.
        Can be input as a uuid or a Protein object.
    :param initial_smiles: The SMILES for the ligand.
    :param num_trajectories: The number of trajectories to run.
    :param simulation_time_ns: How long to run the simulation for, in nanoseconds.
    :param ligand_residue_name: The name of the residue corresponding to the ligand.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """

    if isinstance(protein, Protein):
        protein = protein.uuid

    workflow = stjames.PoseAnalysisMolecularDynamicsWorkflow(
        protein_uuid=protein,
        initial_smiles=initial_smiles,
        num_trajectories=num_trajectories,
        simulation_time_ns=simulation_time_ns,
        ligand_residue_name=ligand_residue_name,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "pose_analysis_md",
        "workflow_data": workflow.model_dump(serialize_as_any=True, mode="json"),
        "initial_smiles": initial_smiles,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


def submit_batch_docking_workflow(
    smiles_list: list[str],
    protein: str | Protein,
    pocket: list[list[float]],
    executable: str = "qvina2",
    scoring_function: str = "vina",
    exhaustiveness: float = 8,
    name: str = "Batch Docking Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a batch docking workflow to the API.

    :param smiles_list: The SMILES strings to dock.
    :param protein: The protein to dock. Can be input as a uuid or a Protein object.
    :param executable: Which docking implementation to use.
    :param scoring_function: Which docking scoring function to use.
    :param exhaustiveness: Which exhaustiveness to employ.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted batch docking workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """

    docking_settings = {
        "executable": executable,
        "exhaustiveness": exhaustiveness,
        "scoring_function": scoring_function,
    }

    workflow = stjames.BatchDockingWorkflow(
        initial_smiles_list=smiles_list,
        target=protein,
        pocket=pocket,
        docking_settings=docking_settings,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "batch_docking",
        "workflow_data": workflow.model_dump(serialize_as_any=True, mode="json"),
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


def submit_msa_workflow(
    initial_protein_sequences: list[stjames.ProteinSequence | str],
    output_formats: list[stjames.MSAFormat],
    name: str = "MSA Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a multiple sequence alignment (MSA) workflow to the API.

    :param initial_protein_sequences: List of protein sequences to align, as ProteinSequence objects
    sor strings.
    :param output_formats: List of output formats for the resulting MSA files.
    :param name: The name to assign to the workflow.
    :param folder_uuid: UUID of the folder where the workflow will be stored.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted MSA workflow.
    :raises HTTPError: If the API request fails.
    """
    workflow = stjames.MSAWorkflow(
        initial_protein_sequences=initial_protein_sequences,
        output_formats=output_formats,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "msa",
        "workflow_data": workflow.model_dump(serialize_as_any=True, mode="json"),
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


def submit_admet_workflow(
    initial_smiles: str,
    name: str = "ADMET Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits an ADMET workflow to the API.

    :param initial_smiles: The molecule used in the workflow.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to store the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """

    workflow = stjames.ADMETWorkflow(initial_smiles=initial_smiles)

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "admet",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_smiles": initial_smiles,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


def submit_membrane_permeability_workflow(
    initial_molecule: dict[str, Any] | StJamesMolecule | RdkitMol | str,
    method: Literal["gnn-mtl", "pypermm"] = "gnn-mtl",
    name: str = "ADMET Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a membrane permeability workflow to the API.

    :param initial_molecule: The molecule used in the workflow.
    :param method: The method used to compute membrane permeability.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to store the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """

    data: dict[str, Any] = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "membrane_permeability",
        "max_credits": max_credits,
    }

    match method:
        case "gnn-mtl":
            assert isinstance(initial_molecule, str)
            workflow = stjames.MembranePermeabilityWorkflow(
                initial_smiles=initial_molecule,
                membrane_permeability_method="chemprop_ohlsson2025",
            )

            data["initial_smiles"] = initial_molecule

        case "pypermm":
            if isinstance(initial_molecule, str):
                raise ValueError("Cannot specify molecule as SMILES for PyPermm")
            elif isinstance(initial_molecule, StJamesMolecule):
                initial_molecule = initial_molecule.model_dump(mode="json")
            elif isinstance(initial_molecule, RdkitMol):
                initial_molecule = StJamesMolecule.from_rdkit(initial_molecule, cid=0).model_dump(
                    mode="json"
                )

            workflow = stjames.MembranePermeabilityWorkflow(
                initial_molecule=initial_molecule,
                membrane_permeability_method="pypermm",
            )

            data["initial_molecule"] = initial_molecule

        case _:
            raise ValueError(f"Unexpected {method=}")

    data["workflow_data"] = workflow.model_dump(mode="json")

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())
