"""Multistage optimization workflow - optimize molecules with staged methods."""

import stjames

from ..calculation import Calculation, retrieve_calculation
from ..molecule import Molecule
from ..types import MoleculeInput, SolventInput
from ..utils import api_client
from .base import (
    Mode,
    Workflow,
    WorkflowResult,
    molecule_to_dict,
    register_result,
)


@register_result("multistage_opt")
class MultiStageOptResult(WorkflowResult):
    """Result from a multistage-optimization workflow."""

    _stjames_class = stjames.MultiStageOptWorkflow

    def __post_init__(self) -> None:
        """Parse workflow data and eagerly fetch final calculation."""
        super().__post_init__()
        uuids = self.calculation_uuids
        if uuids:
            self._cache["final_calculation"] = retrieve_calculation(uuids[-1])

    def __repr__(self) -> str:
        energy = self.energy
        e_str = f"{energy:.6f}" if energy is not None else "None"
        return f"<MultiStageOptResult energy={e_str}>"

    @property
    def calculation_uuids(self) -> list[str]:
        """UUIDs of all calculations in the optimization stages."""
        calcs = getattr(self._workflow, "calculations", None) or []
        return [c for c in calcs if c is not None]

    @property
    def calculations(self) -> list[Calculation]:
        """
        All optimization stage calculations (lazily fetched).

        Typically includes xTB pre-optimization, DFT optimization, and final
        single-point. Access `final_calculation` for just the last stage.
        """
        if "calculations" not in self._cache:
            uuids = self.calculation_uuids
            calcs = []
            for i, uuid in enumerate(uuids):
                # Reuse eagerly-fetched final calculation
                if i == len(uuids) - 1 and "final_calculation" in self._cache:
                    calcs.append(self._cache["final_calculation"])
                else:
                    calcs.append(retrieve_calculation(uuid))
            self._cache["calculations"] = calcs
        return self._cache["calculations"]

    @property
    def final_calculation(self) -> Calculation | None:
        """The final Calculation object with full molecule data (eagerly fetched)."""
        return self._cache.get("final_calculation")

    @property
    def molecule(self) -> Molecule | None:
        """The final optimized molecule geometry."""
        calc = self.final_calculation
        return calc.molecule if calc else None

    @property
    def energy(self) -> float | None:
        """Energy of the final optimized molecule (Hartree)."""
        mol = self.molecule
        return mol.energy if mol else None

    @property
    def charges(self) -> list[float] | None:
        """Partial charges on each atom."""
        mol = self.molecule
        return mol.charges if mol else None

    @property
    def dipole(self) -> tuple[float, float, float] | None:
        """Dipole moment vector (Debye)."""
        mol = self.molecule
        return mol.dipole if mol else None


def submit_multistage_optimization_workflow(
    initial_molecule: MoleculeInput,
    mode: Mode = Mode.RAPID,
    solvent: SolventInput = None,
    xtb_preopt: bool = True,
    transition_state: bool = False,
    frequencies: bool = False,
    name: str = "Multistage Optimization Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a multistage-optimization workflow to the API.

    :param initial_molecule: Molecule to optimize.
    :param mode: Mode to run the calculation in.
    :param solvent: Solvent for the final single-point calculation.
    :param xtb_preopt: Whether to pre-optimize with xTB.
    :param transition_state: Whether this is a transition state optimization.
    :param frequencies: Whether to calculate frequencies.
    :param name: Name of the workflow.
    :param folder_uuid: UUID of the folder to place the workflow in.
    :param max_credits: Maximum number of credits to use for the workflow.
    :returns: Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    mol_dict = molecule_to_dict(initial_molecule)

    workflow = stjames.MultiStageOptWorkflow(
        initial_molecule=mol_dict,
        mode=mode,
        solvent=solvent,
        xtb_preopt=xtb_preopt,
        transition_state=transition_state,
        frequencies=frequencies,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "multistage_opt",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": mol_dict,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())
