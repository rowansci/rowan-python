"""Basic calculation workflow - perform quantum chemical calculations."""

import stjames
from rdkit import Chem

from ..calculation import Calculation, retrieve_calculation
from ..molecule import Molecule
from ..utils import api_client
from .base import (
    MoleculeInput,
    RowanMolecule,
    StJamesMolecule,
    Workflow,
    WorkflowResult,
    register_result,
)
from .constants import HARTREE_TO_KCAL


@register_result("basic_calculation")
class BasicCalculationResult(WorkflowResult):
    """Result from a basic calculation workflow."""

    _stjames_class = stjames.BasicCalculationWorkflow

    def __post_init__(self) -> None:
        """Parse workflow data and eagerly fetch calculation."""
        super().__post_init__()
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
        """The Calculation object with full molecule data."""
        return self._cache.get("calculation")

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

    def get_optimization_energies(self, relative: bool = False) -> list[float | None]:
        """
        Get energies for each optimization step.

        :param relative: If True, return relative energies in kcal/mol (relative to
            the lowest energy step). If False (default), return absolute energies
            in Hartree.
        :return: List of energies for each optimization step.
        """
        energies: list[float | None] = [m.energy for m in self.molecules]

        if relative:
            valid = [e for e in energies if e is not None]
            if valid:
                min_e = min(valid)
                energies = [
                    (e - min_e) * HARTREE_TO_KCAL if e is not None else None for e in energies
                ]

        return energies

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

    if isinstance(initial_molecule, RowanMolecule):
        initial_molecule = initial_molecule._to_stjames().model_dump(mode="json")
    elif isinstance(initial_molecule, StJamesMolecule):
        initial_molecule = initial_molecule.model_dump(mode="json")
    elif isinstance(initial_molecule, Chem.rdchem.Mol | Chem.rdchem.RWMol):
        initial_molecule = stjames.Molecule.from_rdkit(initial_molecule, cid=0)

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
