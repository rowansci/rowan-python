from __future__ import annotations

import cctk
import httpx
from typing import Optional
import numpy as np
import stjames
from dataclasses import dataclass, field
import time

import rowan

API_URL = "https://api.rowansci.com/calculation"

# cf. /opt/miniconda3/envs/openai/lib/python3.9/site-packages/openai/__init__.py

@dataclass
class Client:
    blocking: bool = True
    print: bool = True
    ping_interval: int = 5

    headers: dict = field(init=False)

    def __post_init__(self):
        self.headers = {"X-API-Key": rowan.utils.get_api_key()}

    def compute(
        self,
        input_molecule: cctk.Molecule,
        name: Optional[str] = None,
        folder_id: Optional[str] = None,
        **options,
    ) -> stjames.Calculation | int:
        atomic_numbers = input_molecule.atomic_numbers.view(np.ndarray)
        geometry = input_molecule.geometry.view(np.ndarray)

        atoms = list()
        for i in range(input_molecule.num_atoms()):
            atoms.append(
                stjames.Atom(atomic_number=atomic_numbers[i], position=geometry[i])
            )

        molecule = stjames.Molecule(
            atoms=atoms,
            charge=input_molecule.charge,
            multiplicity=input_molecule.multiplicity,
        )
        settings = stjames.Settings(**options)
        calc = stjames.Calculation(molecules=[molecule], name=name, settings=settings)

        with httpx.Client() as client:
            response = client.post(
                f"{API_URL}",
                headers=self.headers,
                json={
                    "json_data": calc.model_dump(mode="json"),
                    "folder_id": folder_id,
                },
            )
            response.raise_for_status()
            response_dict = response.json()
            calc_id = response_dict["id"]

            if self.print:
                print(f"Calculation {calc_id} submitted")

            if self.blocking:
                while not response_dict["is_finished"]:
                    time.sleep(self.ping_interval)
                    response = client.get(f"{API_URL}/{calc_id}", headers=self.headers)
                    response_dict = response.json()
                    if not response.status_code == 200:
                        raise ValueError(
                            f"Error {response.status_code}: {response.reason}"
                        )

                # we finished! get proper schema
                response = client.get(
                    f"{API_URL}/{calc_id}/stjames", headers=self.headers
                )
                response_dict = response.json()

                if self.print:
                    print(
                        f"Calculation {calc_id} completed after {response_dict['elapsed']:.1f} s of CPU time"
                    )
                    print(
                        f"(View the results at labs.rowansci.com/calculations/{calc_id})"
                    )

                return stjames.Calculation.model_validate(response_dict)

            else:
                return calc_id

    def get(self, calc_id: int) -> stjames.Calculation:
        with httpx.Client() as client:
            response = client.get(f"{API_URL}/{calc_id}/stjames", headers=self.headers)
            response.raise_for_status()
            response_dict = response.json()
            return stjames.Calculation.model_validate(response_dict)

    def stop(self, calc_id: int) -> None:
        with httpx.Client() as client:
            response = client.post(f"{API_URL}/{calc_id}/stop", headers=self.headers)
            response.raise_for_status()
