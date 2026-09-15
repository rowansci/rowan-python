"""Electronic properties workflow - calculate electronic properties."""

from dataclasses import dataclass
from typing import Any

import httpx
import stjames
from pydantic import ValidationError
from stjames import Engine, Method

from ..folder import Folder
from ..utils import api_client
from .base import (
    StructureInput,
    Workflow,
    WorkflowResult,
    molecule_to_dict,
    register_result,
    require_coordinates,
)


@dataclass(frozen=True, slots=True)
class MolecularOrbital:
    """Molecular orbital with cube data.

    Attributes:
        points: cube data as (x, y, z, value) points (Bohr)
        occupation: occupation number (0, 1, or 2)
        energy: orbital energy (Hartree)
    """

    points: tuple[tuple[float, float, float, float], ...]
    occupation: int
    energy: float


@register_result("electronic_properties")
class ElectronicPropertiesResult(WorkflowResult):
    """Result from an electronic-properties workflow."""

    _stjames_class = stjames.ElectronicPropertiesWorkflow

    def __post_init__(self) -> None:
        # Computed results (dipole, charges, cubes, orbitals, ...) live only in
        # the gzipped S3 blob, not in the workflow's object_data row. Fetch and
        # parse it once on construction so accessors below see populated data.
        super().__post_init__()
        if not self.complete:
            return
        try:
            with api_client() as client:
                response = client.get(f"/orbitals/{self.workflow_uuid}/compressed_json")
                response.raise_for_status()
                populated = stjames.ElectronicPropertiesWorkflow.model_validate(response.json())
            object.__setattr__(self, "_workflow", populated)
        except (httpx.HTTPError, ValidationError):
            pass

    def __repr__(self) -> str:
        return f"<ElectronicPropertiesResult dipole={self.dipole} D>"

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
        if not (cube := self._workflow.density_cube):
            return None
        return [(p.x, p.y, p.z, p.val) for p in cube.cube_points]

    @property
    def electrostatic_potential_cube(self) -> list[tuple[float, float, float, float]] | None:
        """Electrostatic potential cube as list of (x, y, z, value) points."""
        if not (cube := self._workflow.electrostatic_potential_cube):
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
    initial_molecule: StructureInput,
    method: Method | str = "b97_3c",
    basis_set: str | None = None,
    engine: Engine | str | None = None,
    compute_density_cube: bool = True,
    compute_electrostatic_potential_cube: bool = True,
    compute_num_occupied_orbitals: int = 1,
    compute_num_virtual_orbitals: int = 1,
    name: str = "Electronic Properties Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
    is_draft: bool = False,
) -> Workflow:
    """Submits an electronic-properties workflow to the API.

    Args:
        initial_molecule: molecule to calculate electronic properties for
        method: method to use for the calculation
        basis_set: basis set to use (if any)
        engine: compute engine, see `Engine`. Auto-selected from method if not specified
        compute_density_cube: whether to compute the density cube
        compute_electrostatic_potential_cube: whether to compute the electrostatic
            potential cube
        compute_num_occupied_orbitals: number of occupied orbitals to save
        compute_num_virtual_orbitals: number of virtual orbitals to save
        name: name of the workflow
        folder_uuid: UUID of the folder to place the workflow in
        folder: destination folder
        max_credits: maximum credits for the workflow
        webhook_url: URL that Rowan will POST to when the workflow completes
        is_draft: save as a draft without starting execution

    Returns:
        submitted workflow

    Raises:
        ValueError: method is not supported by the engine
        httpx.HTTPStatusError: request to the API fails
    """
    require_coordinates(initial_molecule)
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder:
        folder_uuid = folder.uuid
    mol_dict = molecule_to_dict(initial_molecule)

    if isinstance(method, str):
        method = Method(method)
    if isinstance(engine, str):
        engine = Engine(engine)

    # stjames validates method/engine compatibility and auto-selects the engine if None
    settings_kwargs: dict[str, Any] = {"method": method, "basis_set": basis_set}
    if engine is not None:
        settings_kwargs["engine"] = engine
    settings = stjames.Settings(**settings_kwargs)

    workflow = stjames.ElectronicPropertiesWorkflow(
        initial_molecule=mol_dict,
        settings=settings,
        compute_density_cube=compute_density_cube,
        compute_electrostatic_potential_cube=compute_electrostatic_potential_cube,
        compute_num_occupied_orbitals=compute_num_occupied_orbitals,
        compute_num_virtual_orbitals=compute_num_virtual_orbitals,
    )

    data = {
        "workflow_type": "electronic_properties",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": mol_dict,
        "name": name,
        "folder_uuid": folder_uuid,
        "max_credits": max_credits,
        "webhook_url": webhook_url,
        "is_draft": is_draft,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())
