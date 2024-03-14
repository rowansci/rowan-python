from __future__ import annotations

import cctk
import httpx
from typing import Optional
import stjames
from dataclasses import dataclass, field
import time

import rowan

API_URL = "https://api.rowansci.com"


@dataclass
class Client:
    blocking: bool = True
    print: bool = True
    ping_interval: int = 5
    delete_when_finished: bool = False

    headers: dict = field(init=False)

    def __post_init__(self):
        self.headers = {"X-API-Key": rowan.utils.get_api_key()}

    def compute(
        self,
        type: Optional[str] = "calculation",
        input_mol: Optional[cctk.Molecule] = None,
        input_smiles: Optional[str] = None,
        name: Optional[str] = None,
        folder_id: Optional[str] = None,
        **options,
    ) -> dict | str:
        if (input_mol is None) == (input_smiles is None):
            raise ValueError("Must specify exactly one of ``input_smiles`` and ``input_molecule``!")

        if input_mol is None:
            input_mol = cctk.Molecule.new_from_smiles(input_smiles)

        molecule = rowan.utils.cctk_to_stjames(input_mol)

        with httpx.Client() as client:
            if type == "calculation":
                settings = stjames.Settings(**options)
                calc = stjames.Calculation(molecules=[molecule], name=name, settings=settings)

                response = client.post(
                    f"{API_URL}/calculation",
                    headers=self.headers,
                    json={
                        "json_data": calc.model_dump(mode="json"),
                        "folder_id": folder_id,
                    },
                )

            elif type == "pka":
                response = client.post(
                    f"{API_URL}/pka",
                    headers=self.headers,
                    json={
                        "initial_molecule": molecule.model_dump(mode="json"),
                        "name": name,
                        "folder_id": folder_id,
                        "workflow_data": options,
                    },
                )

            response.raise_for_status()
            response_dict = response.json()
            calc_id = response_dict["uuid"]

        if self.blocking:
            while not self.is_finished(calc_id, type):
                time.sleep(self.ping_interval)
            result = self.get(calc_id, type)

            if self.delete_when_finished:
                self.delete(calc_id, type)

            return result

        else:
            return calc_id

    def is_finished(self, calc_id: str, type: str = "calculation") -> bool:
        with httpx.Client() as client:
            if type == "calculation":
                response = client.get(f"{API_URL}/calculation/{calc_id}", headers=self.headers)
                response.raise_for_status()
                response_dict = response.json()
                status = response_dict["status"]

            elif type == "pka":
                response = client.get(f"{API_URL}/pka/{calc_id}", headers=self.headers)
                response.raise_for_status()
                response_dict = response.json()
                status = response_dict["object_status"]

            else:
                raise ValueError(f"Unknown type ``{type}``!")

        return status in [2, 3, 4]

    def get(self, calc_id: str, type: str = "calculation") -> dict:
        with httpx.Client() as client:
            if type == "calculation":
                response = client.get(f"{API_URL}/calculation/{calc_id}/stjames", headers=self.headers)
                response.raise_for_status()
                response_dict = response.json()
                return response_dict

            elif type == "pka":
                response = client.get(f"{API_URL}/pka/{calc_id}", headers=self.headers)
                response.raise_for_status()
                response_dict = response.json()
                return response_dict["object_data"]

            else:
                raise ValueError(f"Unknown type ``{type}``!")

    def stop(self, calc_id: str, type: str = "calculation") -> None:
        with httpx.Client() as client:
            if type == "calculation":
                response = client.post(f"{API_URL}/calculation/{calc_id}/stop", headers=self.headers)
                response.raise_for_status()

            elif type == "pka":
                response = client.post(f"{API_URL}/pka/{calc_id}/stop", headers=self.headers)
                response.raise_for_status()

            else:
                raise ValueError(f"Unknown type ``{type}``!")


    def delete(self, calc_id: str, type: str = "calculation") -> None:
        with httpx.Client() as client:
            if type == "calculation":
                response = client.delete(f"{API_URL}/calculation/{calc_id}", headers=self.headers)
                response.raise_for_status()

            elif type == "pka":
                response = client.delete(f"{API_URL}/folder/{calc_id}", headers=self.headers)
                response.raise_for_status()

            else:
                raise ValueError(f"Unknown type ``{type}``!")

