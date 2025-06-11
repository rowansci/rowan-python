from typing import Any, Optional

import stjames

from .utils import api_client


class Workflow:
    @classmethod
    def submit(
        cls,
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
            return response.json()

    @classmethod
    def retrieve(cls, uuid: stjames.UUID) -> dict[str, Any]:
        with api_client() as client:
            response = client.get(f"/workflow/{uuid}")
            response.raise_for_status()
            return response.json()

    @classmethod
    def update(
        cls,
        uuid: stjames.UUID,
        name: Optional[str] = None,
        parent_uuid: Optional[stjames.UUID] = None,
        notes: Optional[str] = None,
        starred: Optional[bool] = None,
        email_when_complete: Optional[bool] = None,
        public: Optional[bool] = None,
    ) -> None:
        old_data = cls.retrieve(uuid)

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
            return response.json()

    @classmethod
    def status(cls, uuid: stjames.UUID) -> int:
        return cls.retrieve(uuid)["object_status"]

    @classmethod
    def is_finished(cls, uuid: stjames.UUID) -> bool:
        status = cls.status(uuid)

        return status in {
            stjames.Status.COMPLETED_OK.value,
            stjames.Status.FAILED.value,
            stjames.Status.STOPPED.value,
        }

    @classmethod
    def stop(cls, uuid: stjames.UUID) -> None:
        with api_client() as client:
            response = client.post(f"/workflow/{uuid}/stop")
            response.raise_for_status()

    @classmethod
    def delete(cls, uuid: stjames.UUID) -> None:
        with api_client() as client:
            response = client.delete(f"/workflow/{uuid}")
            response.raise_for_status()

    @classmethod
    def delete_data(cls, uuid: stjames.UUID) -> None:
        with api_client() as client:
            response = client.delete(f"/workflow/{uuid}/delete_workflow_data")
            response.raise_for_status()

    @classmethod
    def list(
        cls,
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
            return response.json()
