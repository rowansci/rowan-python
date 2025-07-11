import time
from datetime import datetime
from typing import Any, Self, TypeAlias

import stjames
from pydantic import BaseModel, Field
from rdkit import Chem

from .utils import api_client

RdkitMol: TypeAlias = Chem.rdchem.Mol | Chem.rdchem.RWMol


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
    data: dict[str, Any] = Field(default={}, alias="object_data")
    email_when_complete: bool
    max_credits: int | None = None
    elapsed: float | None = None
    credits_charged: float

    class Config:  # noqa: D106
        validate_by_name = True

    def __repr__(self) -> str:
        return f"<Workflow name='{self.name}' created_at='{self.created_at}'>"

    def fetch_latest(self, in_place: bool = False) -> Self:
        """
        Loads workflow data from the database and updates the current instance.

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


def submit_workflow(
    workflow_type: str,
    workflow_data: dict[str, Any] | None = None,
    initial_molecule: dict[str, Any] | stjames.Molecule | RdkitMol | None = None,
    initial_smiles: str | None = None,
    name: str | None = None,
    folder_uuid: str | None = None,
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
    :return: A Workflow object representing the submitted workflow.
    :raises ValueError: If neither `initial_smiles` nor a valid `initial_molecule` is provided.
    :raises HTTPError: If the API request fails.
    """
    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": workflow_type,
        "workflow_data": workflow_data,
    }

    if initial_smiles is not None:
        data["initial_smiles"] = initial_smiles
    elif isinstance(initial_molecule, stjames.Molecule):
        data["initial_molecule"] = initial_molecule.model_dump()
    elif isinstance(initial_molecule, dict):
        data["initial_molecule"] = initial_molecule
    elif isinstance(initial_molecule, RdkitMol):
        data["initial_molecule"] = stjames.Molecule.from_rdkit(initial_molecule, cid=0).model_dump()
    else:
        raise ValueError("You must provide either `initial_smiles` or a valid `initial_molecule`.")

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


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


def retrieve_calculation_molecules(uuid: str) -> list[dict[str, Any]]:
    """
    Retrieves a list of molecules from a calculation.

    :param uuid: The UUID of the calculation to retrieve molecules from.
    :return: A list of dictionaries representing the molecules in the calculation.
    :raises HTTPError: If the API request fails.
    """
    with api_client() as client:
        response = client.get(f"/calculation/{uuid}/molecules")
        response.raise_for_status()
        return response.json()


def list_workflows(
    parent_uuid: str | None = None,
    name_contains: str | None = None,
    public: bool | None = None,
    starred: bool | None = None,
    status: int | None = None,
    workflow_type: str | None = None,
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
    initial_molecule: dict[str, Any] | stjames.Molecule | RdkitMol,
    method: stjames.Method | str = "uma_m_omol",
    tasks: list[str] | None = None,
    mode: str = "auto",
    engine: str = "omol25",
    name: str = "Basic Calculation Workflow",
    folder_uuid: str | None = None,
) -> Workflow:
    """
    Submit a basic calculation workflow to the API.

    :param initial_molecule: The molecule to perform the calculation on.
    :param method: The method to use for the calculation.
    See [list of available methods](https://github.com/rowansci/stjames-public/blob/master/stjames/method.py)
    for options.
    :param tasks: A list of tasks to perform for the calculation.
    :param mode: The mode to run the calculation in. See [list of available modes](https://github.com/rowansci/stjames-public/blob/master/stjames/mode.py)
    for options.
    :param engine: The engine to use for the calculation. See [list of available engines](https://github.com/rowansci/stjames-public/blob/master/stjames/engine.py)
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if not tasks:
        tasks = ["optimize"]

    if isinstance(initial_molecule, stjames.Molecule):
        initial_molecule = initial_molecule.model_dump()
    elif isinstance(initial_molecule, RdkitMol):
        initial_molecule = stjames.Molecule.from_rdkit(initial_molecule, cid=0)

    if isinstance(method, str):
        method = stjames.Method(method)

    workflow_data = {
        "settings": {
            "method": method.name,
            "tasks": tasks,
            "mode": mode,
        },
        "engine": engine,
    }

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "basic_calculation",
        "workflow_data": workflow_data,
        "initial_molecule": initial_molecule,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


def submit_conformer_search_workflow(
    initial_molecule: dict[str, Any] | stjames.Molecule | RdkitMol,
    conf_gen_mode: str = "rapid",
    final_method: stjames.Method | str = "aimnet2_wb97md3",
    solvent: str | None = None,
    transistion_state: bool = False,
    name: str = "Conformer Search Workflow",
    folder_uuid: str | None = None,
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
    :param transistion_state: Whether to optimize the transition state.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if isinstance(initial_molecule, stjames.Molecule):
        initial_molecule = initial_molecule.model_dump()
    elif isinstance(initial_molecule, RdkitMol):
        initial_molecule = stjames.Molecule.from_rdkit(initial_molecule, cid=0)

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
        opt_settings={"transition_state": transistion_state, "constraints": []},
    )

    msos = stjames.MultiStageOptSettings(
        mode=stjames.Mode.MANUAL,
        xtb_preopt=True,
        optimization_settings=[opt_settings],
    )

    workflow_data = {
        "multistage_opt_settings": msos.model_dump(),
        "conf_gen_mode": conf_gen_mode,
        "mso_mode": "manual",
        "solvent": solvent,
        "transistion_state": transistion_state,
    }

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "conformer_search",
        "workflow_data": workflow_data,
        "initial_molecule": initial_molecule,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


def submit_pka_workflow(
    initial_molecule: dict[str, Any] | stjames.Molecule | RdkitMol,
    pka_range: tuple[int, int] = (2, 12),
    deprotonate_elements: list[int] | None = None,
    protonate_elements: list[int] | None = None,
    mode: str = "careful",
    name: str = "pKa Workflow",
    folder_uuid: str | None = None,
) -> Workflow:
    """
    Submits a pKa workflow to the API.

    :param initial_molecule: The molecule to calculate the pKa of.
    :param pka_range: The range of pKa values to calculate.
    :param deprotonate_elements: The elements to deprotonate. Given by atomic number.
    :param protonate_elements: The elements to protonate. Given by atomic number.
    :param mode: The mode to run the calculation in. See
    [list of available modes](https://github.com/rowansci/stjames-public/blob/master/stjames/mode.py)
    for options.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if isinstance(initial_molecule, stjames.Molecule):
        initial_molecule = initial_molecule.model_dump()
    elif isinstance(initial_molecule, RdkitMol):
        initial_molecule = stjames.Molecule.from_rdkit(initial_molecule, cid=0)

    protonate_elements = protonate_elements or [7]
    deprotonate_elements = deprotonate_elements or [7, 8, 16]

    workflow_data = {
        "pka_range": pka_range,
        "deprotonate_elements": deprotonate_elements,
        "protonate_elements": protonate_elements,
        "mode": mode,
    }

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "pka",
        "workflow_data": workflow_data,
        "initial_molecule": initial_molecule,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


def submit_scan_workflow(
    initial_molecule: dict[str, Any] | stjames.Molecule | RdkitMol,
    scan_settings: stjames.ScanSettings | dict[str, Any] | None = None,
    calculation_engine: str = "omol25",
    calculation_method: stjames.Method | str = "uma_m_omol",
    wavefront_propagation: bool = True,
    name: str = "Scan Workflow",
    folder_uuid: str | None = None,
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
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if isinstance(initial_molecule, stjames.Molecule):
        initial_molecule = initial_molecule.model_dump()
    elif isinstance(initial_molecule, RdkitMol):
        initial_molecule = stjames.Molecule.from_rdkit(initial_molecule, cid=0)

    if isinstance(calculation_method, str):
        calculation_method = stjames.Method(calculation_method)

    workflow_data = {
        "wavefront_propagation": wavefront_propagation,
        "scan_settings": scan_settings,
        "calc_engine": calculation_engine,
        "calc_settings": {
            "method": calculation_method.name,
            "corrections": [],
            "tasks": ["optimize"],
            "mode": "auto",
            "opt_settings": {"constraints": []},
        },
    }

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "scan",
        "workflow_data": workflow_data,
        "initial_molecule": initial_molecule,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


def submit_irc_workflow(
    initial_molecule: dict[str, Any] | stjames.Molecule | RdkitMol | None = None,
    method: stjames.Method | str = "uma_m_omol",
    engine: str = "omol25",
    preopt: bool = True,
    step_size: float = 0.05,
    max_irc_steps: int = 30,
    name: str = "IRC Workflow",
    folder_uuid: str | None = None,
) -> Workflow:
    """
    Submits an Intrinsic Reaction Coordinate (IRC) workflow to the API.

    :param initial_molecule: The initial molecule to perform the IRC calculation on.
    :param method: The computational method to use for the IRC calculation.
    See [list of available methods](https://github.com/rowansci/stjames-public/blob/master/stjames/method.py)
    for options.
    :param engine: The computational engine to use for the calculation. See [list of available engines](https://github.com/rowansci/stjames-public/blob/master/stjames/engine.py)
    :param preopt: Whether to perform a pre-optimization of the molecule.
    :param step_size: The step size to use for the IRC calculation.
    :param max_irc_steps: The maximum number of IRC steps to perform.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :return: A Workflow object representing the submitted IRC workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """

    if isinstance(initial_molecule, stjames.Molecule):
        initial_molecule = initial_molecule.model_dump()
    elif isinstance(initial_molecule, RdkitMol):
        initial_molecule = stjames.Molecule.from_rdkit(initial_molecule, cid=0)

    if isinstance(method, str):
        method = stjames.Method(method)

    workflow_data = {
        "settings": {
            "method": method.name,
            "tasks": [],
            "corrections": [],
            "mode": "auto",
        },
        "engine": engine,
        "preopt": preopt,
        "step_size": step_size,
        "max_irc_steps": max_irc_steps,
        "mode": "manual",
    }

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "irc",
        "workflow_data": workflow_data,
        "initial_molecule": initial_molecule,
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
    name: str = "Cofolding Workflow",
    model: str = stjames.CofoldingModel.BOLTZ_2.value,
    folder_uuid: str | None = None,
) -> Workflow:
    """
    Submits a protein cofolding workflow to the API.

    :param initial_protein_sequences: The sequences of the proteins to be cofolded.
    :param initial_smiles_list: A list of SMILES strings for the ligands to be cofolded with.
    :param ligand_binding_affinity_index: The index of the ligand for which to compute the binding affinity.
    :param use_msa_server: Whether to use the MSA server for the computation.
    :param use_potentials: Whether to use potentials for the computation.
    :param name: The name of the workflow.
    :param model: The model to use for the computation.
    :param folder_uuid: The UUID of the folder to store the workflow in.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """  # noqa: E501
    workflow_data = {
        "use_msa_server": use_msa_server,
        "use_potentials": use_potentials,
        "model": model,
        "ligand_binding_affinity_index": ligand_binding_affinity_index,
        "initial_smiles_list": initial_smiles_list,
        "initial_protein_sequences": initial_protein_sequences,
    }
    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "protein_cofolding",
        "workflow_data": workflow_data,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())
