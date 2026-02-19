"""BDE workflow - Bond Dissociation Energy calculations."""

from dataclasses import dataclass
from typing import Any

import stjames
from rdkit import Chem

from ..utils import api_client
from .base import RdkitMol, StJamesMolecule, Workflow, WorkflowResult, register_result


@dataclass(frozen=True, slots=True)
class BDEEntry:
    """A bond dissociation energy result."""

    fragment_idxs: tuple[int, ...]
    energy: float | None = None


@register_result("bde")
class BDEResult(WorkflowResult):
    """Result from a Bond-Dissociation Energy (BDE) workflow."""

    _stjames_class = stjames.BDEWorkflow

    def __repr__(self) -> str:
        n = len(self.bdes)
        return f"<BDEResult energy={self.energy} bdes={n}>"

    @property
    def energy(self) -> float | None:
        """Energy of the molecule (Hartree)."""
        return self._workflow.optimization_energy

    @property
    def bdes(self) -> list[BDEEntry]:
        """Bond dissociation energies."""
        return [
            BDEEntry(
                fragment_idxs=tuple(b.fragment_idxs),
                energy=b.energy,
            )
            for b in self._workflow.bdes
        ]


def submit_bde_workflow(
    initial_molecule: dict[str, Any] | StJamesMolecule | RdkitMol,
    mode: str = "rapid",
    atoms: list[int] | None = None,
    all_CH: bool = False,
    all_CX: bool = False,
    name: str = "BDE Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a bond dissociation energy (BDE) workflow to the API.

    :param initial_molecule: The molecule to calculate BDEs for.
    :param mode: The mode to run the calculation in.
    :param atoms: List of atom indices (1-indexed) to dissociate.
    :param all_CH: Whether to dissociate all C-H bonds.
    :param all_CX: Whether to dissociate all C-X bonds (X = halogen).
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if isinstance(initial_molecule, StJamesMolecule):
        initial_molecule = initial_molecule.model_dump(mode="json")
    elif isinstance(initial_molecule, Chem.rdchem.Mol | Chem.rdchem.RWMol):
        initial_molecule = StJamesMolecule.from_rdkit(initial_molecule, cid=0).model_dump(
            mode="json"
        )

    workflow = stjames.BDEWorkflow(
        initial_molecule=initial_molecule,
        mode=mode,
        atoms=atoms or [],
        all_CH=all_CH,
        all_CX=all_CX,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "bde",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": initial_molecule,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


__all__ = ["BDEEntry", "BDEResult", "submit_bde_workflow"]
