"""Calculation results from Rowan workflows."""

from typing import Self

import stjames
from pydantic import BaseModel, ConfigDict

from .utils import api_client


class Calculation(BaseModel):
    """
    A Rowan calculation result.

    :ivar uuid: The UUID of the calculation
    :ivar name: The name of the calculation
    :ivar status: The status code of the calculation
    :ivar elapsed: Execution time in seconds
    :ivar engine: The compute engine used
    :ivar molecules: List of molecules (geometries) from this calculation
    """

    uuid: str
    name: str | None = None
    status: int | None = None
    elapsed: float | None = None
    engine: str | None = None
    molecules: list[stjames.Molecule] = []

    model_config = ConfigDict(arbitrary_types_allowed=True)

    def __repr__(self) -> str:
        n_mols = len(self.molecules)
        energy = self.energy
        e_str = f"{energy:.6f}" if energy is not None else "None"
        return f"<Calculation energy={e_str} molecules={n_mols} uuid='{self.uuid}'>"

    @property
    def energy(self) -> float | None:
        """Energy of the final molecule (Hartree)."""
        if self.molecules:
            return self.molecules[-1].energy
        return None

    @property
    def molecule(self) -> stjames.Molecule | None:
        """The final molecule geometry."""
        if self.molecules:
            return self.molecules[-1]
        return None

    def refresh(self, in_place: bool = True) -> Self:
        """
        Fetch the latest calculation data from the API.

        :param in_place: If True, update this instance. If False, return new instance.
        :return: The updated Calculation object.
        """
        with api_client() as client:
            response = client.get(f"/calculation/{self.uuid}/stjames")
            response.raise_for_status()
            data = response.json()

        if not in_place:
            return self.__class__(uuid=self.uuid, **data)

        self.name = data.get("name")
        self.status = data.get("status")
        self.elapsed = data.get("elapsed")
        self.engine = data.get("engine")
        self.molecules = [stjames.Molecule.model_validate(m) for m in data.get("molecules", [])]
        return self


def retrieve_calculation(uuid: str) -> Calculation:
    """
    Retrieve a calculation from the API by UUID.

    :param uuid: The UUID of the calculation.
    :return: A Calculation object with the fetched data.
    :raises requests.HTTPError: If the API request fails.
    """
    with api_client() as client:
        response = client.get(f"/calculation/{uuid}/stjames")
        response.raise_for_status()
        data = response.json()

    molecules = [stjames.Molecule.model_validate(m) for m in data.get("molecules", [])]

    return Calculation(
        uuid=uuid,
        name=data.get("name"),
        status=data.get("status"),
        elapsed=data.get("elapsed"),
        engine=data.get("engine"),
        molecules=molecules,
    )


__all__ = ["Calculation", "retrieve_calculation"]
