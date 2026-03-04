"""Tautomer search workflow - find tautomeric forms of molecules."""

from dataclasses import dataclass
from typing import Any

import stjames
from rdkit import Chem

from ..calculation import retrieve_calculation
from ..molecule import Molecule
from ..utils import api_client
from .base import (
    Message,
    Mode,
    RdkitMol,
    StJamesMolecule,
    Workflow,
    WorkflowResult,
    parse_messages,
    register_result,
)


@dataclass(frozen=True, slots=True)
class Tautomer:
    """A tautomer result."""

    energy: float
    """Energy in Hartree."""

    weight: float
    """Boltzmann weight (sum to 1.0 across all tautomers)."""

    predicted_relative_energy: float
    """Relative energy in kcal/mol (relative to lowest energy tautomer)."""

    structure_uuids: tuple[str, ...]
    """UUIDs of the structure calculations."""


@register_result("tautomers")
class TautomerResult(WorkflowResult):
    """Result from a tautomer search workflow."""

    _stjames_class = stjames.TautomerWorkflow

    def __repr__(self) -> str:
        tautomers = self.tautomers
        n = len(tautomers)
        if tautomers:
            lowest = min(tautomers, key=lambda t: t.energy)
            return f"<TautomerResult tautomers={n} lowest_energy={lowest.energy:.6f}>"
        return f"<TautomerResult tautomers={n}>"

    @property
    def tautomers(self) -> list[Tautomer]:
        """List of tautomers with energies and weights."""
        return [
            Tautomer(
                energy=t.energy,
                weight=t.weight,
                predicted_relative_energy=t.predicted_relative_energy,
                structure_uuids=tuple(s.uuid for s in t.structures if s.uuid),
            )
            for t in self._workflow.tautomers
        ]

    @property
    def messages(self) -> list[Message]:
        """Any messages or warnings from the workflow."""
        return parse_messages(getattr(self._workflow, "messages", None))

    @property
    def molecules(self) -> list[Molecule]:
        """Molecules for all tautomers.

        Note: Makes one API call per tautomer on first access.
        Results are cached. Call clear_cache() to refresh.
        """
        if "all_molecules" not in self._cache:
            mols = []
            for t in self.tautomers:
                if t.structure_uuids:
                    calc = retrieve_calculation(t.structure_uuids[0])
                    if calc.molecule:
                        mols.append(calc.molecule)
            self._cache["all_molecules"] = mols
        return self._cache["all_molecules"]


def submit_tautomer_search_workflow(
    initial_molecule: dict[str, Any] | stjames.Molecule | RdkitMol,
    mode: Mode = Mode.CAREFUL,
    name: str = "Tautomer Search Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a tautomer search workflow to the API.

    :param initial_molecule: The molecule to find tautomers for.
    :param mode: The mode to run the calculation in (reckless, rapid, careful).
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: If the request to the API fails.
    """
    if isinstance(initial_molecule, StJamesMolecule):
        initial_molecule = initial_molecule.model_dump(mode="json")
    elif isinstance(initial_molecule, Chem.rdchem.Mol | Chem.rdchem.RWMol):
        initial_molecule = stjames.Molecule.from_rdkit(initial_molecule, cid=0)

    workflow = stjames.TautomerWorkflow(
        initial_molecule=initial_molecule,
        mode=mode,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "tautomers",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": initial_molecule,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


__all__ = ["Tautomer", "TautomerResult", "submit_tautomer_search_workflow"]
