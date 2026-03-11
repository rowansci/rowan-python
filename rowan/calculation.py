"""Calculation results from Rowan workflows."""

from typing import Self

import stjames
from pydantic import BaseModel, ConfigDict

from .molecule import Molecule
from .utils import api_client


def _parse_molecules(data: list[dict]) -> list[Molecule]:
    """Parse molecule dicts into Molecule objects."""
    return [Molecule._from_stjames(stjames.Molecule.model_validate(m)) for m in data]


class Calculation(BaseModel):
    """
    Rowan calculation result.

    :param uuid: UUID of the calculation.
    :param name: Name of the calculation.
    :param status: Status code of the calculation.
    :param elapsed: Execution time in seconds.
    :param engine: Compute engine used.
    :param molecules: Molecules (geometries) from this calculation.
    """

    uuid: str
    name: str | None = None
    status: int | None = None
    elapsed: float | None = None
    engine: str | None = None
    _molecules: list[Molecule] = []

    model_config = ConfigDict(arbitrary_types_allowed=True)

    def __init__(self, molecules: list[Molecule] | None = None, **kwargs):
        super().__init__(**kwargs)
        object.__setattr__(self, "_molecules", molecules or [])

    def __repr__(self) -> str:
        n_mols = len(self._molecules)
        energy = self.energy
        e_str = f"{energy:.6f}" if energy is not None else "None"
        return f"<Calculation energy={e_str} molecules={n_mols} uuid='{self.uuid}'>"

    @property
    def molecules(self) -> list[Molecule]:
        """List of molecules (geometries) from this calculation."""
        return self._molecules

    @property
    def energy(self) -> float | None:
        """Energy of the final molecule (Hartree)."""
        if self._molecules:
            return self._molecules[-1].energy
        return None

    @property
    def molecule(self) -> Molecule | None:
        """The final molecule geometry."""
        if self._molecules:
            return self._molecules[-1]
        return None

    def refresh(self, in_place: bool = True) -> Self:
        """
        Fetch the latest calculation data from the API.

        :param in_place: If True, update this instance in-place. If False, return new instance.
        :returns: Updated Calculation object.
        """
        with api_client() as client:
            response = client.get(f"/calculation/{self.uuid}/stjames")
            response.raise_for_status()
            data = response.json()

        molecules = _parse_molecules(data.get("molecules", []))

        if not in_place:
            return self.__class__(
                uuid=self.uuid,
                name=data.get("name"),
                status=data.get("status"),
                elapsed=data.get("elapsed"),
                engine=data.get("engine"),
                molecules=molecules,
            )

        self.name = data.get("name")
        self.status = data.get("status")
        self.elapsed = data.get("elapsed")
        self.engine = data.get("engine")
        object.__setattr__(self, "_molecules", molecules)
        return self


def retrieve_calculation(uuid: str) -> Calculation:
    """
    Retrieve a calculation from the API by UUID.

    :param uuid: UUID of the calculation to retrieve.
    :returns: Calculation object with the fetched data.
    :raises requests.HTTPError: If the API request fails.
    """
    with api_client() as client:
        response = client.get(f"/calculation/{uuid}/stjames")
        response.raise_for_status()
        data = response.json()

    molecules = _parse_molecules(data.get("molecules", []))

    return Calculation(
        uuid=uuid,
        name=data.get("name"),
        status=data.get("status"),
        elapsed=data.get("elapsed"),
        engine=data.get("engine"),
        molecules=molecules,
    )


__all__ = ["Calculation", "retrieve_calculation"]
