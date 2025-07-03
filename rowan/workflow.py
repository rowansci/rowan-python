from datetime import datetime
from typing import Any, Self, TypeAlias

import stjames
from pydantic import BaseModel, Field
from rdkit import Chem

from .utils import api_client

RdkitMol: TypeAlias = Chem.rdchem.Mol | Chem.rdchem.RWMol


class Workflow(BaseModel):
    name: str
    uuid: str
    created_at: datetime | None = None
    updated_at: datetime | None = None
    started_at: datetime | None = None
    completed_at: datetime | None = None
    status: int | None = Field(default=None, alias="object_status")
    parent_uuid: str
    notes: str | None = None
    starred: bool | None = None
    public: bool | None = None
    workflow_type: str | None = Field(default=None, alias="object_type")
    data: dict[str, Any] | None = Field(default=None, alias="object_data")
    email_when_complete: bool | None = None
    max_credits: int | None = None
    elapsed: float | None = None
    credits_charged: float | None = None

    class Config:
        """
        Configuration for pydantic
        """

        validate_by_name = True

    def __repr__(self) -> str:
        return f"<Workflow name='{self.name}' created_at='{self.created_at}'>"

    def load_data(self) -> Self:
        with api_client() as client:
            response = client.get(f"/workflow/{self.uuid}")
            response.raise_for_status()
            return type(self)(**response.json())

    def update(
        self,
        name: str | None = None,
        parent_uuid: str | None = None,
        notes: str | None = None,
        starred: bool | None = None,
        email_when_complete: bool | None = None,
        public: bool | None = None,
    ) -> None:
        old_data = self.load_data()

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
            updated = response.json()

        # Update self with new data
        for key, value in updated.items():
            setattr(self, key, value)

    def get_status(self) -> stjames.Status:
        return stjames.Status(self.load_data().status or 0)

    def is_finished(self) -> bool:
        status = self.get_status()

        return status in {
            stjames.Status.COMPLETED_OK,
            stjames.Status.FAILED,
            stjames.Status.STOPPED,
        }

    def stop(self) -> None:
        with api_client() as client:
            response = client.post(f"/workflow/{self.uuid}/stop")
            response.raise_for_status()

    def delete(self) -> None:
        with api_client() as client:
            response = client.delete(f"/workflow/{self.uuid}")
            response.raise_for_status()

    def delete_data(self) -> None:
        with api_client() as client:
            response = client.delete(f"/workflow/{self.uuid}/delete_workflow_data")
            response.raise_for_status()


def submit_workflow(
    workflow_type: str,
    workflow_data: dict[str, Any],
    initial_molecule: dict[str, Any] | stjames.Molecule | RdkitMol | None = None,
    initial_smiles: str | None = None,
    name: str | None = None,
    folder_uuid: str | None = None,
) -> Workflow:
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
    with api_client() as client:
        response = client.get(f"/workflow/{uuid}")
        response.raise_for_status()
        return Workflow(**response.json())


def retrieve_calculation_molecules(uuid: str) -> list[dict[str, Any]]:
    with api_client() as client:
        response = client.get(f"/calculation/{uuid}/molecules")
        response.raise_for_status()
        return response.json()


def update(
    uuid: str,
    name: str | None = None,
    parent_uuid: str | None = None,
    notes: str | None = None,
    starred: bool | None = None,
    email_when_complete: bool | None = None,
    public: bool | None = None,
) -> Workflow:
    new_data: dict[str, Any] = {}

    if name is not None:
        new_data["name"] = name

    if parent_uuid is not None:
        new_data["parent_uuid"] = parent_uuid

    if notes is not None:
        new_data["notes"] = notes

    if starred is not None:
        new_data["starred"] = starred

    if email_when_complete is not None:
        new_data["email_when_complete"] = email_when_complete

    if public is not None:
        new_data["public"] = public

    with api_client() as client:
        response = client.post(f"/workflow/{uuid}", json=new_data)
        response.raise_for_status()
        return Workflow(**response.json())


def get_status(uuid: str) -> stjames.Status:
    return stjames.Status(retrieve_workflow(uuid).status or 0)


def is_finished(uuid: str) -> bool:
    status = get_status(uuid)

    return status in {
        stjames.Status.COMPLETED_OK,
        stjames.Status.FAILED,
        stjames.Status.STOPPED,
    }


def stop(uuid: str) -> None:
    with api_client() as client:
        response = client.post(f"/workflow/{uuid}/stop")
        response.raise_for_status()


def delete(uuid: str) -> None:
    with api_client() as client:
        response = client.delete(f"/workflow/{uuid}")
        response.raise_for_status()


def delete_data(uuid: str) -> None:
    with api_client() as client:
        response = client.delete(f"/workflow/{uuid}/delete_workflow_data")
        response.raise_for_status()


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
    initial_molecule: dict[str, Any] | stjames.Molecule | RdkitMol | None = None,
    method: stjames.Method | str = "uma_m_omol",
    tasks: list[str] | None = None,
    mode: str = "auto",
    engine: str = "omol25",
    name: str = "Basic Calculation Workflow",
    folder_uuid: str | None = None,
) -> Workflow:
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
    initial_molecule: dict[str, Any] | stjames.Molecule | RdkitMol | None = None,
    conf_gen_mode: str = "rapid",
    final_method: stjames.Method | str = "aimnet2_wb97md3",
    solvent: str | None = None,
    transistion_state: bool = False,
    name: str = "Conformer Search Workflow",
    folder_uuid: str | None = None,
) -> Workflow:
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
    initial_molecule: dict[str, Any] | stjames.Molecule | RdkitMol | None = None,
    pka_range: tuple[int, int] = (2, 12),
    deprotonate_elements: list[int] | None = None,
    protonate_elements: list[int] | None = None,
    mode: str = "careful",
    name: str = "pKa Workflow",
    folder_uuid: str | None = None,
) -> Workflow:
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
    initial_molecule: dict[str, Any] | stjames.Molecule | RdkitMol | None = None,
    scan_settings: stjames.ScanSettings | dict[str, Any] | None = None,
    calculation_engine: str = "omol25",
    calculation_method: stjames.Method | str = "uma_m_omol",
    wavefront_propagation: bool = True,
    name: str = "Scan Workflow",
    folder_uuid: str | None = None,
) -> Workflow:
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
