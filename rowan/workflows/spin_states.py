"""Spin-states workflow - calculate energies of different spin multiplicities."""

from dataclasses import dataclass
from typing import Any

import stjames

from ..calculation import Calculation, retrieve_calculation
from ..folder import Folder
from ..types import MoleculeInput, SolventInput
from ..utils import api_client
from .base import (
    Message,
    Mode,
    Workflow,
    WorkflowResult,
    molecule_to_dict,
    parse_messages,
    register_result,
)
from .constants import to_relative_kcal


def _validate_multiplicity(mol_dict: dict[str, Any], multiplicity: int) -> None:
    """
    Validate that a spin multiplicity is compatible with the molecule.

    Uses stjames.Molecule.check_electron_sanity() for validation.

    :param mol_dict: Molecule dict with atomic_numbers and charge.
    :param multiplicity: Spin multiplicity to validate.
    :raises ValueError: If multiplicity is invalid for this molecule.
    """
    # Create a copy with the test multiplicity and validate using stjames
    test_dict = {**mol_dict, "multiplicity": multiplicity}
    mol = stjames.Molecule(**test_dict)
    mol.check_electron_sanity()  # type: ignore[operator]


@dataclass(frozen=True, slots=True)
class SpinState:
    """
    Spin state result.

    :param multiplicity: Spin multiplicity (1=singlet, 2=doublet, 3=triplet, etc.).
    :param energy: Energy in Hartree.
    :param calculation_uuids: UUIDs for each optimization stage (for multistage optimization).
    """

    multiplicity: int
    energy: float
    calculation_uuids: tuple[str | None, ...]


@register_result("spin_states")
class SpinStatesResult(WorkflowResult):
    """Result from a spin-states workflow."""

    _stjames_class = stjames.SpinStatesWorkflow

    def __repr__(self) -> str:
        states = self.spin_states
        n = len(states)
        if states:
            g = min(states, key=lambda s: s.energy)
            return f"<SpinStatesResult states={n} ground=(mult={g.multiplicity}, E={g.energy} H)>"
        return f"<SpinStatesResult states={n}>"

    @property
    def spin_states(self) -> list[SpinState]:
        """List of spin states with energies, in submission order."""
        return [
            SpinState(
                multiplicity=ss.multiplicity,
                energy=ss.energy,
                calculation_uuids=tuple(ss.calculation) if ss.calculation else (),
            )
            for ss in self._workflow.spin_states
        ]

    @property
    def messages(self) -> list[Message]:
        """Any messages or warnings from the workflow."""
        return parse_messages(getattr(self._workflow, "messages", None))

    def get_calculation(self, multiplicity: int, stage: int = -1) -> Calculation:
        """
        Fetch the calculation for a specific spin state.

        .. note::
            Makes one API call per spin state on first access.
            Results are cached. Call clear_cache() to refresh.

        :param multiplicity: Spin multiplicity to fetch.
        :param stage: Optimization stage (-1 for final stage).
        :returns: Calculation object with molecule and energy data.
        :raises ValueError: If the multiplicity is not found or has no calculation.
        """
        for state in self.spin_states:
            if state.multiplicity == multiplicity:
                if not state.calculation_uuids:
                    raise ValueError(f"Spin state {multiplicity} has no calculations")
                uuid = state.calculation_uuids[stage]
                if uuid is None:
                    raise ValueError(
                        f"Spin state {multiplicity} has no calculation at stage {stage}"
                    )

                cache_key = f"spin_state_{multiplicity}_{stage}"
                if cache_key not in self._cache:
                    self._cache[cache_key] = retrieve_calculation(uuid)
                return self._cache[cache_key]

        raise ValueError(f"Spin state with multiplicity {multiplicity} not found")

    def get_energies(self, relative: bool = False) -> list[float]:
        """
        Get energies for each spin state.

        :param relative: If True, return relative energies in kcal/mol (relative to
            the ground state / lowest energy spin state). If False (default),
            return absolute energies in Hartree.
        :returns: List of energies for each spin state.
        """
        energies: list[float] = [s.energy for s in self.spin_states]
        return to_relative_kcal(energies) if relative else energies


def submit_spin_states_workflow(
    initial_molecule: MoleculeInput,
    states: list[int],
    mode: Mode = Mode.RAPID,
    solvent: SolventInput = None,
    xtb_preopt: bool = True,
    frequencies: bool = False,
    name: str = "Spin States Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
) -> Workflow:
    """
    Submits a spin-states workflow to the API.

    :param initial_molecule: Molecule to calculate spin states for.
    :param states: List of multiplicities to calculate
        (e.g., [1, 3, 5] for singlet, triplet, quintet).
    :param mode: Mode to run the calculation in.
    :param solvent: Solvent to use for the calculation.
    :param xtb_preopt: Whether to pre-optimize with xTB.
    :param frequencies: Whether to calculate frequencies.
    :param name: Name of the workflow.
    :param folder_uuid: UUID of the folder to place the workflow in.
    :param folder: Folder object to store the workflow in.
    :param max_credits: Maximum number of credits to use for the workflow.
    :param webhook_url: URL that Rowan will POST to when the workflow completes.
    :returns: Workflow object representing the submitted workflow.
    :raises ValueError: If any multiplicity is incompatible with the molecule.
    :raises requests.HTTPError: If the request to the API fails.
    """
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder:
        folder_uuid = folder.uuid
    mol_dict = molecule_to_dict(initial_molecule)

    # Validate that all multiplicities are compatible with the molecule
    for mult in states:
        _validate_multiplicity(mol_dict, mult)

    workflow = stjames.SpinStatesWorkflow(
        initial_molecule=mol_dict,
        states=states,
        mode=mode,
        solvent=solvent,
        xtb_preopt=xtb_preopt,
        frequencies=frequencies,
    )

    data = {
        "workflow_type": "spin_states",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": mol_dict,
        "name": name,
        "folder_uuid": folder_uuid,
        "max_credits": max_credits,
        "webhook_url": webhook_url,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())
