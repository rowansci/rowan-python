"""Basic calculation workflow - perform quantum chemical calculations."""

import stjames

from ..calculation import Calculation, retrieve_calculation
from ..folder import Folder
from ..molecule import Molecule
from ..utils import api_client
from .base import (
    MoleculeInput,
    Workflow,
    WorkflowResult,
    molecule_to_dict,
    register_result,
)
from .constants import to_relative_kcal


@register_result("basic_calculation")
class BasicCalculationResult(WorkflowResult):
    """Result from a basic-calculation workflow."""

    _stjames_class = stjames.BasicCalculationWorkflow

    def __post_init__(self) -> None:
        super().__post_init__()
        if self.eager:
            calc_uuid = getattr(self._workflow, "calculation_uuid", None)
            if calc_uuid:
                self._cache["calculation"] = retrieve_calculation(calc_uuid)

    def __repr__(self) -> str:
        energy = self.energy
        e_str = f"{energy:.6f}" if energy is not None else "None"
        return f"<BasicCalculationResult energy={e_str}>"

    @property
    def calculation_uuid(self) -> str | None:
        """The UUID of the calculation."""
        return getattr(self._workflow, "calculation_uuid", None)

    @property
    def calculation(self) -> Calculation | None:
        """The Calculation object with full molecule data (lazily fetched)."""
        if "calculation" not in self._cache:
            calc_uuid = self.calculation_uuid
            if calc_uuid:
                self._cache["calculation"] = retrieve_calculation(calc_uuid)
        return self._cache.get("calculation")

    @property
    def molecule(self) -> Molecule | None:
        """The final molecule geometry with all computed properties."""
        calc = self.calculation
        return calc.molecule if calc else None

    @property
    def molecules(self) -> list[Molecule]:
        """All molecules from the calculation (e.g., optimization trajectory)."""
        calc = self.calculation
        return calc.molecules if calc else []

    @property
    def energy(self) -> float | None:
        """Energy of the final molecule (Hartree)."""
        mol = self.molecule
        return mol.energy if mol else None

    def optimization_energies(self, relative: bool = False) -> list[float]:
        """
        Energies for each optimization step.

        :param relative: If True, return relative energies in kcal/mol (relative to
            the lowest energy step). If False (default), return absolute energies
            in Hartree.
        :returns: List of energies for each optimization step.
        """
        energies: list[float] = [m.energy for m in self.molecules if m.energy is not None]
        return to_relative_kcal(energies) if relative else energies

    @property
    def charges(self) -> list[float] | None:
        """Partial charges on each atom."""
        mol = self.molecule
        return mol.charges if mol else None

    @property
    def spin_densities(self) -> list[float] | None:
        """Spin densities on each atom (for open-shell systems)."""
        mol = self.molecule
        return mol.spin_densities if mol else None

    @property
    def dipole(self) -> tuple[float, float, float] | None:
        """Dipole moment vector (Debye)."""
        mol = self.molecule
        return mol.dipole if mol else None

    @property
    def frequencies(self) -> list[float] | None:
        """Vibrational frequencies in cm^-1 (if frequency calculation was performed)."""
        mol = self.molecule
        return mol.frequencies if mol else None


def submit_basic_calculation_workflow(
    initial_molecule: MoleculeInput,
    method: stjames.Method | str = "omol25_conserving_s",
    basis_set: stjames.BasisSet | str | None = None,
    tasks: list[str] | None = None,
    mode: str = "auto",
    engine: str | None = None,
    name: str = "Basic Calculation Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submit a basic-calculation workflow to the API.

    :param initial_molecule: Molecule to perform the calculation on.
    :param method: Method to use for the calculation.
    :param basis_set: Basis set to use (if any).
    :param tasks: Tasks to perform for the calculation.
    :param mode: Mode to run the calculation in.
    :param engine: Engine to use for the calculation.
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
    if not tasks:
        tasks = ["optimize"]

    initial_molecule = molecule_to_dict(initial_molecule)

    if isinstance(method, str):
        method = stjames.Method(method)

    workflow = stjames.BasicCalculationWorkflow(
        initial_molecule=initial_molecule,
        settings=stjames.Settings(
            method=method,
            basis_set=basis_set,
            tasks=tasks,
            mode=mode,
        ),
        engine=engine or method.default_engine(),
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "basic_calculation",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": initial_molecule,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())
