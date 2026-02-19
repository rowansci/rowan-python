"""Basic calculation workflow - perform quantum chemical calculations."""

from typing import Any

import stjames
from rdkit import Chem

from ..calculation import Calculation, retrieve_calculation
from ..molecule import Molecule
from ..utils import api_client
from .base import RdkitMol, StJamesMolecule, Workflow, WorkflowResult, register_result


@register_result("basic_calculation")
class BasicCalculationResult(WorkflowResult):
    """Result from a basic calculation workflow."""

    _stjames_class = stjames.BasicCalculationWorkflow

    def __post_init__(self) -> None:
        """Parse workflow data and eagerly fetch calculation."""
        super().__post_init__()
        # Eagerly fetch calculation data so molecule is immediately available
        calc_uuid = getattr(self._workflow, "calculation_uuid", None)
        if calc_uuid:
            self._cache["calculation"] = retrieve_calculation(calc_uuid)

    def __repr__(self) -> str:
        energy = self.energy
        e_str = f"{energy:.6f}" if energy is not None else "None"
        parts = [f"energy={e_str}"]
        if self.solvent:
            parts.append(f"solvent={self.solvent!r}")
        return f"<BasicCalculationResult {' '.join(parts)}>"

    @property
    def calculation_uuid(self) -> str | None:
        """UUID of the calculation."""
        return getattr(self._workflow, "calculation_uuid", None)

    @property
    def calculation(self) -> Calculation | None:
        """The Calculation object with full molecule data."""
        return self._cache.get("calculation")

    @property
    def method(self) -> stjames.Method | None:
        """The computational method used."""
        settings = getattr(self._workflow, "settings", None)
        return settings.method if settings else None

    @property
    def engine(self) -> str | None:
        """The compute engine used."""
        return getattr(self._workflow, "engine", None)

    @property
    def solvent(self) -> str | None:
        """The solvent used (if any)."""
        settings = getattr(self._workflow, "settings", None)
        if settings and settings.solvent_settings:
            return settings.solvent_settings.get("solvent")
        return None

    @property
    def tasks(self) -> list[str] | None:
        """The tasks performed (e.g., ['optimize'], ['energy'], ['frequency'])."""
        settings = getattr(self._workflow, "settings", None)
        return list(settings.tasks) if settings and settings.tasks else None

    @property
    def molecule(self) -> Molecule | None:
        """The final molecule geometry with all computed properties."""
        calc = self.calculation
        return calc.molecule if calc else None

    @property
    def molecules(self) -> list:
        """All molecules from the calculation (e.g., optimization trajectory)."""
        calc = self.calculation
        return calc.molecules if calc else []

    @property
    def energy(self) -> float | None:
        """Energy of the final molecule (Hartree)."""
        mol = self.molecule
        return mol.energy if mol else None

    @property
    def optimization_energies(self) -> list[float]:
        """Energies for each optimization step (if optimization was performed)."""
        return [m.energy for m in self.molecules if m.energy is not None]

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
    initial_molecule: dict[str, Any] | StJamesMolecule | RdkitMol,
    method: stjames.Method | str = "uma_m_omol",
    basis_set: stjames.BasisSet | str | None = None,
    tasks: list[str] | None = None,
    mode: str = "auto",
    engine: str | None = None,
    name: str = "Basic Calculation Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submit a basic calculation workflow to the API.

    :param initial_molecule: The molecule to perform the calculation on.
    :param method: The method to use for the calculation.
    :param basis_set: The basis_set to use (if any).
    :param tasks: A list of tasks to perform for the calculation.
    :param mode: The mode to run the calculation in.
    :param engine: The engine to use for the calculation.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if not tasks:
        tasks = ["optimize"]

    if isinstance(initial_molecule, StJamesMolecule):
        initial_molecule = initial_molecule.model_dump(mode="json")
    elif isinstance(initial_molecule, Chem.rdchem.Mol | Chem.rdchem.RWMol):
        initial_molecule = StJamesMolecule.from_rdkit(initial_molecule, cid=0)

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


__all__ = ["BasicCalculationResult", "submit_basic_calculation_workflow"]
