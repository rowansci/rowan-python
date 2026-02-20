"""Tautomer search workflow - find tautomeric forms of molecules."""

from dataclasses import dataclass
from typing import Any

import stjames
from rdkit import Chem

from ..utils import api_client
from .base import RdkitMol, StJamesMolecule, Workflow, WorkflowResult, register_result


@dataclass(frozen=True, slots=True)
class Tautomer:
    """A tautomer result."""

    energy: float
    weight: float | None = None
    predicted_relative_energy: float | None = None


@register_result("tautomers")
class TautomerResult(WorkflowResult):
    """Result from a tautomer search workflow."""

    _stjames_class = stjames.TautomerWorkflow

    def __repr__(self) -> str:
        tautomers = self.tautomers
        n = len(tautomers)
        if tautomers:
            lowest = min(tautomers, key=lambda t: t.energy)
            return f"<TautomerResult tautomers={n} lowest_energy={lowest.energy}>"
        return f"<TautomerResult tautomers={n}>"

    @property
    def tautomers(self) -> list[Tautomer]:
        """List of tautomer structures."""
        return [
            Tautomer(
                energy=t.energy,
                weight=t.weight,
                predicted_relative_energy=t.predicted_relative_energy,
            )
            for t in self._workflow.tautomers
        ]


def submit_tautomer_search_workflow(
    initial_molecule: dict[str, Any] | stjames.Molecule | RdkitMol,
    mode: str = "careful",
    name: str = "Tautomer Search Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a tautomer search workflow to the API.

    :param initial_molecule: The molecule to calculate the tautomers of.
    :param mode: The mode to run the calculation in.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
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
