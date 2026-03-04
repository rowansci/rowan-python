"""Electronic properties workflow - calculate electronic properties."""

from dataclasses import dataclass
from typing import Any

import stjames
from rdkit import Chem

from ..utils import api_client
from .base import RdkitMol, StJamesMolecule, Workflow, WorkflowResult, register_result


@dataclass(frozen=True, slots=True)
class MolecularOrbital:
    """A molecular orbital with cube data."""

    points: tuple[tuple[float, float, float, float], ...]  # (x, y, z, value)
    occupation: int
    energy: float


@register_result("electronic_properties")
class ElectronicPropertiesResult(WorkflowResult):
    """Result from an electronic properties workflow."""

    _stjames_class = stjames.ElectronicPropertiesWorkflow

    def __repr__(self) -> str:
        return f"<ElectronicPropertiesResult dipole={self.dipole}>"

    @property
    def dipole(self) -> tuple[float, float, float] | None:
        """Dipole moment vector (Debye)."""
        return self._workflow.dipole

    @property
    def quadrupole(
        self,
    ) -> (
        tuple[
            tuple[float, float, float],
            tuple[float, float, float],
            tuple[float, float, float],
        ]
        | None
    ):
        """Quadrupole moment tensor (Debye·Å)."""
        return self._workflow.quadrupole

    @property
    def mulliken_charges(self) -> list[float] | None:
        """Mulliken partial charges on each atom."""
        return self._workflow.mulliken_charges

    @property
    def lowdin_charges(self) -> list[float] | None:
        """Löwdin partial charges on each atom."""
        return self._workflow.lowdin_charges

    @property
    def wiberg_bond_orders(self) -> list[tuple[int, int, float]]:
        """Wiberg bond orders as (atom1, atom2, order) tuples."""
        return self._workflow.wiberg_bond_orders

    @property
    def mayer_bond_orders(self) -> list[tuple[int, int, float]]:
        """Mayer bond orders as (atom1, atom2, order) tuples."""
        return self._workflow.mayer_bond_orders

    @property
    def density_cube(self) -> list[tuple[float, float, float, float]] | None:
        """Electron density cube as list of (x, y, z, value) points."""
        cube = self._workflow.density_cube
        if cube is None:
            return None
        return [(p.x, p.y, p.z, p.val) for p in cube.cube_points]

    @property
    def electrostatic_potential_cube(self) -> list[tuple[float, float, float, float]] | None:
        """Electrostatic potential cube as list of (x, y, z, value) points."""
        cube = self._workflow.electrostatic_potential_cube
        if cube is None:
            return None
        return [(p.x, p.y, p.z, p.val) for p in cube.cube_points]

    @property
    def molecular_orbitals(self) -> dict[int, MolecularOrbital]:
        """Molecular orbitals indexed by orbital number."""
        return {
            k: MolecularOrbital(
                points=tuple((p.x, p.y, p.z, p.val) for p in v.cube_points),
                occupation=v.occupation,
                energy=v.energy,
            )
            for k, v in self._workflow.molecular_orbitals.items()
        }


def submit_electronic_properties_workflow(
    initial_molecule: dict[str, Any] | stjames.Molecule | RdkitMol,
    method: stjames.Method | str = "b97_3c",
    basis_set: str | None = None,
    compute_density_cube: bool = True,
    compute_electrostatic_potential_cube: bool = True,
    compute_num_occupied_orbitals: int = 1,
    compute_num_virtual_orbitals: int = 1,
    name: str = "Electronic Properties Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits an electronic properties workflow to the API.

    :param initial_molecule: The molecule to calculate electronic properties for.
    :param method: The method to use for the calculation.
    :param basis_set: The basis set to use (if any).
    :param compute_density_cube: Whether to compute the density cube.
    :param compute_electrostatic_potential_cube: Whether to compute the electrostatic
        potential cube.
    :param compute_num_occupied_orbitals: Number of occupied orbitals to save.
    :param compute_num_virtual_orbitals: Number of virtual orbitals to save.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if isinstance(initial_molecule, StJamesMolecule):
        initial_molecule = initial_molecule.model_dump(mode="json")
    elif isinstance(initial_molecule, Chem.rdchem.Mol | Chem.rdchem.RWMol):
        initial_molecule = stjames.Molecule.from_rdkit(initial_molecule, cid=0).model_dump(
            mode="json"
        )

    if isinstance(method, str):
        method = stjames.Method(method)

    settings = stjames.Settings(method=method, basis_set=basis_set)

    workflow = stjames.ElectronicPropertiesWorkflow(
        initial_molecule=initial_molecule,
        settings=settings,
        compute_density_cube=compute_density_cube,
        compute_electrostatic_potential_cube=compute_electrostatic_potential_cube,
        compute_num_occupied_orbitals=compute_num_occupied_orbitals,
        compute_num_virtual_orbitals=compute_num_virtual_orbitals,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "electronic_properties",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": initial_molecule,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


__all__ = [
    "ElectronicPropertiesResult",
    "MolecularOrbital",
    "submit_electronic_properties_workflow",
]
