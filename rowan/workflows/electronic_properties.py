"""Electronic properties workflow - calculate electronic properties."""

from dataclasses import dataclass

import stjames

from ..folder import Folder
from ..types import MoleculeInput
from ..utils import api_client
from .base import Workflow, WorkflowResult, molecule_to_dict, register_result


@dataclass(frozen=True, slots=True)
class MolecularOrbital:
    """Molecular orbital with cube data.

    :param points: Cube data as (x, y, z, value) points (Bohr).
    :param occupation: Occupation number (0, 1, or 2).
    :param energy: Orbital energy (Hartree).
    """

    points: tuple[tuple[float, float, float, float], ...]
    occupation: int
    energy: float


@register_result("electronic_properties")
class ElectronicPropertiesResult(WorkflowResult):
    """Result from an electronic-properties workflow."""

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
        """Quadrupole moment tensor (Debye*A)."""
        return self._workflow.quadrupole

    @property
    def mulliken_charges(self) -> list[float] | None:
        """Mulliken partial charges on each atom."""
        return self._workflow.mulliken_charges

    @property
    def lowdin_charges(self) -> list[float] | None:
        """Lowdin partial charges on each atom."""
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

    @property
    def homo(self) -> MolecularOrbital | None:
        """Highest occupied molecular orbital (HOMO). Energy in Hartree."""
        orbs = self.molecular_orbitals
        occupied = [idx for idx, mo in orbs.items() if mo.occupation > 0]
        return orbs[max(occupied)] if occupied else None

    @property
    def lumo(self) -> MolecularOrbital | None:
        """Lowest unoccupied molecular orbital (LUMO). Energy in Hartree."""
        orbs = self.molecular_orbitals
        virtual = [idx for idx, mo in orbs.items() if mo.occupation == 0]
        return orbs[min(virtual)] if virtual else None

    @property
    def homo_lumo_gap(self) -> float | None:
        """HOMO-LUMO gap (Hartree)."""
        homo = self.homo
        lumo = self.lumo
        if homo is None or lumo is None:
            return None
        return lumo.energy - homo.energy


def submit_electronic_properties_workflow(
    initial_molecule: MoleculeInput,
    method: stjames.Method | str = "b97_3c",
    basis_set: str | None = None,
    compute_density_cube: bool = True,
    compute_electrostatic_potential_cube: bool = True,
    compute_num_occupied_orbitals: int = 1,
    compute_num_virtual_orbitals: int = 1,
    name: str = "Electronic Properties Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits an electronic-properties workflow to the API.

    :param initial_molecule: Molecule to calculate electronic properties for.
    :param method: Method to use for the calculation.
    :param basis_set: Basis set to use (if any).
    :param compute_density_cube: Whether to compute the density cube.
    :param compute_electrostatic_potential_cube: Whether to compute the electrostatic
        potential cube.
    :param compute_num_occupied_orbitals: Number of occupied orbitals to save.
    :param compute_num_virtual_orbitals: Number of virtual orbitals to save.
    :param name: Name of the workflow.
    :param folder_uuid: UUID of the folder to place the workflow in.
    :param folder: Folder object to store the workflow in.
    :param max_credits: Maximum number of credits to use for the workflow.
    :returns: Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if folder is not None and folder_uuid is not None:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder is not None:
        folder_uuid = folder.uuid
    initial_molecule = molecule_to_dict(initial_molecule)

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
