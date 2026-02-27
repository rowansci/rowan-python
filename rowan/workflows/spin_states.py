"""Spin states workflow - calculate energies of different spin multiplicities."""

from dataclasses import dataclass

import stjames

from ..utils import api_client
from .base import (
    Mode,
    MoleculeInput,
    SolventInput,
    Workflow,
    WorkflowResult,
    molecule_to_dict,
    register_result,
)


@dataclass(frozen=True, slots=True)
class SpinState:
    """A spin state result."""

    multiplicity: int
    energy: float


@register_result("spin_states")
class SpinStatesResult(WorkflowResult):
    """Result from a spin states workflow."""

    _stjames_class = stjames.SpinStatesWorkflow

    def __repr__(self) -> str:
        states = self.spin_states
        n = len(states)
        if states:
            g = min(states, key=lambda s: s.energy)
            return f"<SpinStatesResult states={n} ground=(mult={g.multiplicity}, E={g.energy})>"
        return f"<SpinStatesResult states={n}>"

    @property
    def spin_states(self) -> list[SpinState]:
        """List of spin states with energies."""
        return [
            SpinState(
                multiplicity=ss.multiplicity,
                energy=ss.energy,
            )
            for ss in self._workflow.spin_states
        ]

    def get_energies(self, relative: bool = False) -> list[float | None]:
        """
        Get energies for each spin state.

        :param relative: If True, return relative energies in kcal/mol (relative to
            the ground state / lowest energy spin state). If False (default),
            return absolute energies in Hartree.
        :return: List of energies for each spin state.
        """
        energies: list[float | None] = list(self._workflow.energies)

        if relative:
            valid = [e for e in energies if e is not None]
            if valid:
                min_e = min(valid)
                hartree_to_kcal = 627.509
                energies = [
                    (e - min_e) * hartree_to_kcal if e is not None else None
                    for e in energies
                ]

        return energies


def submit_spin_states_workflow(
    initial_molecule: MoleculeInput,
    states: list[int],
    mode: Mode = Mode.RAPID,
    solvent: SolventInput = None,
    xtb_preopt: bool = True,
    frequencies: bool = False,
    name: str = "Spin States Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a spin states workflow to the API.

    :param initial_molecule: The molecule to calculate spin states for.
    :param states: List of multiplicities to calculate
        (e.g., [1, 3, 5] for singlet, triplet, quintet).
    :param mode: The mode to run the calculation in.
    :param solvent: The solvent to use for the calculation.
    :param xtb_preopt: Whether to pre-optimize with xTB.
    :param frequencies: Whether to calculate frequencies.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    mol_dict = molecule_to_dict(initial_molecule)

    workflow = stjames.SpinStatesWorkflow(
        initial_molecule=mol_dict,
        states=states,
        mode=mode,
        solvent=solvent,
        xtb_preopt=xtb_preopt,
        frequencies=frequencies,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "spin_states",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": mol_dict,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


__all__ = ["SpinState", "SpinStatesResult", "submit_spin_states_workflow"]
