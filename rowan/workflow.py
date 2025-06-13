from dataclasses import dataclass
from typing import Any, Optional

import stjames

from .utils import api_client


@dataclass
class Workflow:
    uuid: str
    name: str
    workflow_type: str
    status: int | None

    def retrieve(self) -> "Workflow":
        with api_client() as client:
            response = client.get(f"/workflow/{self.uuid}")
            response.raise_for_status()
            # TODO: return Workflow
            return response.json()

    def update(
        self,
        uuid: stjames.UUID,
        name: Optional[str] = None,
        parent_uuid: Optional[stjames.UUID] = None,
        notes: Optional[str] = None,
        starred: Optional[bool] = None,
        email_when_complete: Optional[bool] = None,
        public: Optional[bool] = None,
    ) -> None:
        old_data = self.retrieve()

        new_data = {}
        new_data["name"] = name if name is not None else old_data["name"]
        new_data["parent_uuid"] = (
            parent_uuid if parent_uuid is not None else old_data["parent_uuid"]
        )
        new_data["notes"] = notes if notes is not None else old_data["notes"]
        new_data["starred"] = starred if starred is not None else old_data["starred"]
        new_data["email_when_complete"] = (
            email_when_complete
            if email_when_complete is not None
            else old_data["email_when_complete"]
        )
        new_data["public"] = public if public is not None else old_data["public"]

        with api_client() as client:
            response = client.post(f"/workflow/{uuid}", json=new_data)
            response.raise_for_status()
            # TODO: return Workflow
            return response.json()

    def status(self) -> int:
        # TODO: fix this access
        return self.retrieve()["object_status"]

    def is_finished(self) -> bool:
        status = self.status()

        return status in {
            stjames.Status.COMPLETED_OK.value,
            stjames.Status.FAILED.value,
            stjames.Status.STOPPED.value,
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
    initial_molecule: dict | stjames.Molecule | None = None,
    initial_smiles: str | None = None,
    name: str | None = None,
    folder_uuid: stjames.UUID | None = None,
) -> dict[str, Any]:
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
    else:
        raise ValueError(
            "You must provide either `initial_smiles` or a valid `initial_molecule`."
        )

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        # TODO: return Workflow
        return response.json()

def retrieve_workflow(uuid: stjames.UUID) -> dict[str, Any]:
    with api_client() as client:
        response = client.get(f"/workflow/{uuid}")
        response.raise_for_status()
        # TODO: return Workflow
        return response.json()

def update(
    uuid: stjames.UUID,
    name: Optional[str] = None,
    parent_uuid: Optional[stjames.UUID] = None,
    notes: Optional[str] = None,
    starred: Optional[bool] = None,
    email_when_complete: Optional[bool] = None,
    public: Optional[bool] = None,
) -> None:
    old_data = retrieve_workflow(uuid)

    new_data = {}
    new_data["name"] = name if name is not None else old_data["name"]
    new_data["parent_uuid"] = (
        parent_uuid if parent_uuid is not None else old_data["parent_uuid"]
    )
    new_data["notes"] = notes if notes is not None else old_data["notes"]
    new_data["starred"] = starred if starred is not None else old_data["starred"]
    new_data["email_when_complete"] = (
        email_when_complete
        if email_when_complete is not None
        else old_data["email_when_complete"]
    )
    new_data["public"] = public if public is not None else old_data["public"]

    with api_client() as client:
        response = client.post(f"/workflow/{uuid}", json=new_data)
        response.raise_for_status()
        # TODO: return Workflow
        return response.json()

def get_status(uuid: stjames.UUID) -> int:
    return retrieve_workflow(uuid)["object_status"]

def is_finished(uuid: stjames.UUID) -> bool:
    status = get_status(uuid)

    return status in {
        stjames.Status.COMPLETED_OK.value,
        stjames.Status.FAILED.value,
        stjames.Status.STOPPED.value,
    }

def stop(uuid: stjames.UUID) -> None:
    with api_client() as client:
        response = client.post(f"/workflow/{uuid}/stop")
        response.raise_for_status()

def delete(uuid: stjames.UUID) -> None:
    with api_client() as client:
        response = client.delete(f"/workflow/{uuid}")
        response.raise_for_status()

def delete_data(uuid: stjames.UUID) -> None:
    with api_client() as client:
        response = client.delete(f"/workflow/{uuid}/delete_workflow_data")
        response.raise_for_status()

def list_workflows(
    parent_uuid: Optional[stjames.UUID] = None,
    name_contains: Optional[str] = None,
    public: Optional[bool] = None,
    starred: Optional[bool] = None,
    object_status: Optional[int] = None,
    object_type: Optional[str] = None,
    page: int = 0,
    size: int = 10,
) -> dict[str, Any]:
    params: dict[str, Any] = {"page": page, "size": size}

    if parent_uuid is not None:
        params["parent_uuid"] = parent_uuid

    if name_contains is not None:
        params["name_contains"] = name_contains

    if public is not None:
        params["public"] = public

    if starred is not None:
        params["starred"] = starred

    if object_status is not None:
        params["object_status"] = object_status

    if object_type is not None:
        params["object_type"] = object_type

    with api_client() as client:
        response = client.get("/workflow", params=params)
        response.raise_for_status()
        # TODO: return Workflows
        return response.json()


def submit_protein_cofolding_workflow(
        initial_protein_sequences: list[str],
        initial_smiles_list: list[str] | None = None,
        ligand_binding_affinity_index: int | None = None,
        use_msa_server: bool = True,
        use_potentials: bool = False,
        name: str = "Cofolding Workflow",
        model: str = stjames.CofoldingModel.BOLTZ_2.value,
        folder_uuid: stjames.UUID | None = None,
    ) -> dict[str, Any]:
        workflow_data = {
            "use_msa_server": use_msa_server,
            "use_potentials": use_potentials,
            "model": model,
            "ligand_binding_affinity_index": ligand_binding_affinity_index,
            "initial_smiles_list": initial_smiles_list,
            "initial_protein_sequences": initial_protein_sequences
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
            # TODO: return Workflow
            return response.json()
