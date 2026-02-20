"""IRC workflow - Intrinsic Reaction Coordinate calculations."""

from typing import Any

import stjames
from rdkit import Chem

from ..calculation import Calculation, retrieve_calculation
from ..molecule import Molecule
from ..utils import api_client
from .base import RdkitMol, StJamesMolecule, Workflow, WorkflowResult, register_result


@register_result("irc")
class IRCResult(WorkflowResult):
    """Result from an Intrinsic Reaction Coordinate (IRC) workflow."""

    _stjames_class = stjames.IRCWorkflow

    def __post_init__(self) -> None:
        """Parse workflow data and eagerly fetch TS calculation."""
        super().__post_init__()
        # Eagerly fetch TS calculation only
        ts_uuid = self._workflow.starting_TS
        if ts_uuid:
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
        """The transition state Calculation (eagerly loaded)."""
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

        :return: List of Calculation objects along the forward path.

        Note: Makes one API call per step. Results are cached.
        """
        if "forward_calcs" not in self._cache:
            calcs = [retrieve_calculation(uuid) for uuid in self.forward_uuids]
            self._cache["forward_calcs"] = calcs
        return self._cache["forward_calcs"]

    def get_backward_calculations(self) -> list["Calculation"]:
        """
        Fetch all backward IRC path calculations.

        :return: List of Calculation objects along the backward path.

        Note: Makes one API call per step. Results are cached.
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

    @property
    def forward_energies(self) -> list[float]:
        """Energies along the forward IRC path (Hartree, lazily fetched)."""
        return [m.energy for m in self.forward_molecules if m.energy is not None]

    @property
    def backward_energies(self) -> list[float]:
        """Energies along the backward IRC path (Hartree, lazily fetched)."""
        return [m.energy for m in self.backward_molecules if m.energy is not None]


def submit_irc_workflow(
    initial_molecule: dict[str, Any] | stjames.Molecule | RdkitMol | None = None,
    method: stjames.Method | str = "uma_m_omol",
    solvent: str | None = None,
    preopt: bool = True,
    step_size: float = 0.05,
    max_irc_steps: int = 30,
    name: str = "IRC Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits an Intrinsic Reaction Coordinate (IRC) workflow to the API.

    :param initial_molecule: The initial molecule to perform the IRC calculation on.
    :param method: The computational method to use for the IRC calculation.
    :param solvent: The solvent to use for the calculation.
    :param preopt: Whether to perform a pre-optimization of the molecule.
    :param step_size: The step size to use for the IRC calculation.
    :param max_irc_steps: The maximum number of IRC steps to perform.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted IRC workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if isinstance(initial_molecule, StJamesMolecule):
        initial_molecule = initial_molecule.model_dump(mode="json")
    elif isinstance(initial_molecule, Chem.rdchem.Mol | Chem.rdchem.RWMol):
        initial_molecule = stjames.Molecule.from_rdkit(initial_molecule, cid=0)

    if isinstance(method, str):
        method = stjames.Method(method)

    workflow = stjames.IRCWorkflow(
        initial_molecule=initial_molecule,
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
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "irc",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": initial_molecule,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


__all__ = ["IRCResult", "submit_irc_workflow"]
