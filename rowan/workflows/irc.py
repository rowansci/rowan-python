"""IRC workflow - Intrinsic Reaction Coordinate calculations."""

import stjames

from ..calculation import Calculation, retrieve_calculation
from ..folder import Folder
from ..molecule import Molecule
from ..types import MoleculeInput, SolventInput
from ..utils import api_client
from .base import (
    Workflow,
    WorkflowResult,
    molecule_to_dict,
    register_result,
)
from .constants import to_relative_kcal


@register_result("irc")
class IRCResult(WorkflowResult):
    """Result from an Intrinsic Reaction Coordinate (IRC) workflow."""

    _stjames_class = stjames.IRCWorkflow

    def __post_init__(self) -> None:
        super().__post_init__()
        if self.complete:
            if ts_uuid := self._workflow.starting_TS:
                self._cache["ts_calc"] = retrieve_calculation(ts_uuid)

    def __repr__(self) -> str:
        n_fwd = len(self.forward_molecules)
        n_bwd = len(self.backward_molecules)
        return f"<IRCResult forward_steps={n_fwd} backward_steps={n_bwd}>"

    @property
    def ts_uuid(self) -> str | None:
        """UUID of the transition state calculation."""
        return self._workflow.starting_TS

    @property
    def forward_uuids(self) -> list[str]:
        """UUIDs of the forward IRC path calculations."""
        return [u for u in (self._workflow.irc_forward or []) if u]

    @property
    def backward_uuids(self) -> list[str]:
        """UUIDs of the backward IRC path calculations."""
        return [u for u in (self._workflow.irc_backward or []) if u]

    @property
    def ts_calculation(self) -> "Calculation | None":
        """The transition state Calculation (lazily fetched)."""
        if "ts_calc" not in self._cache:
            if ts_uuid := self.ts_uuid:
                self._cache["ts_calc"] = retrieve_calculation(ts_uuid)
        return self._cache.get("ts_calc")

    @property
    def ts_molecule(self) -> Molecule | None:
        """The optimized transition state molecule."""
        calc = self.ts_calculation
        return calc.molecule if calc else None

    @property
    def ts_energy(self) -> float | None:
        """Energy of the transition state (Hartree)."""
        mol = self.ts_molecule
        return mol.energy if mol else None

    def get_forward_calculations(self) -> list["Calculation"]:
        """
        Fetch all forward IRC path calculations.

        :returns: List of Calculation objects along the forward path.

        .. note::
            Makes one API call per step. Results are cached.
        """
        if "forward_calcs" not in self._cache:
            calcs = [retrieve_calculation(uuid) for uuid in self.forward_uuids]
            self._cache["forward_calcs"] = calcs
        return self._cache["forward_calcs"]

    def get_backward_calculations(self) -> list["Calculation"]:
        """
        Fetch all backward IRC path calculations.

        :returns: List of Calculation objects along the backward path.

        .. note::
            Makes one API call per step. Results are cached.
        """
        if "backward_calcs" not in self._cache:
            calcs = [retrieve_calculation(uuid) for uuid in self.backward_uuids]
            self._cache["backward_calcs"] = calcs
        return self._cache["backward_calcs"]

    @property
    def forward_molecules(self) -> list[Molecule]:
        """Molecules along the forward IRC path (lazily fetched)."""
        calcs = self.get_forward_calculations()
        return [c.molecule for c in calcs if c.molecule]

    @property
    def backward_molecules(self) -> list[Molecule]:
        """Molecules along the backward IRC path (lazily fetched)."""
        calcs = self.get_backward_calculations()
        return [c.molecule for c in calcs if c.molecule]

    def get_forward_energies(self, relative: bool = False) -> list[float]:
        """
        Get energies along the forward IRC path.

        :param relative: If True, return relative energies in kcal/mol (relative to
            the lowest energy point). If False (default), return absolute energies
            in Hartree.
        :returns: List of energies along the forward path.
        """
        energies: list[float] = [m.energy for m in self.forward_molecules if m.energy is not None]
        return to_relative_kcal(energies) if relative else energies

    def get_backward_energies(self, relative: bool = False) -> list[float]:
        """
        Get energies along the backward IRC path.

        :param relative: If True, return relative energies in kcal/mol (relative to
            the lowest energy point). If False (default), return absolute energies
            in Hartree.
        :returns: List of energies along the backward path.
        """
        energies: list[float] = [m.energy for m in self.backward_molecules if m.energy is not None]
        return to_relative_kcal(energies) if relative else energies


def submit_irc_workflow(
    initial_molecule: MoleculeInput,
    method: stjames.Method | str = "uma_m_omol",
    solvent: SolventInput = None,
    preopt: bool = True,
    step_size: float = 0.05,
    max_irc_steps: int = 30,
    name: str = "IRC Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
) -> Workflow:
    """
    Submits an Intrinsic Reaction Coordinate (IRC) workflow to the API.

    :param initial_molecule: Transition state molecule to start IRC from.
    :param method: Computational method to use for the IRC calculation.
    :param solvent: Solvent to use for the calculation.
    :param preopt: Whether to perform a pre-optimization of the TS.
    :param step_size: Step size to use for the IRC calculation.
    :param max_irc_steps: Maximum number of IRC steps to perform.
    :param name: Name of the workflow.
    :param folder_uuid: UUID of the folder to place the workflow in.
    :param folder: Folder object to store the workflow in.
    :param max_credits: Maximum number of credits to use for the workflow.
    :param webhook_url: URL that Rowan will POST to when the workflow completes.
    :returns: Workflow object representing the submitted IRC workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder:
        folder_uuid = folder.uuid
    mol_dict = molecule_to_dict(initial_molecule)

    if isinstance(method, str):
        method = stjames.Method(method)

    workflow = stjames.IRCWorkflow(
        initial_molecule=mol_dict,
        settings=stjames.Settings(
            method=method,
            tasks=[],
            corrections=[],
            mode="auto",
        ),
        solvent=solvent,
        preopt=preopt,
        step_size=step_size,
        max_irc_steps=max_irc_steps,
        mode="manual",
    )

    data = {
        "workflow_type": "irc",
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
