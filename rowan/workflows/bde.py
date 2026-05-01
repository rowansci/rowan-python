"""BDE workflow - Bond Dissociation Energy calculations."""

from dataclasses import dataclass

import stjames
from stjames.workflows.bde import find_AB_bonds as _find_AB_bonds
from stjames.workflows.bde import find_CH_bonds as _find_CH_bonds
from stjames.workflows.bde import find_CX_bonds as _find_CX_bonds

from ..folder import Folder
from ..utils import api_client
from .base import (
    MoleculeInput,
    Workflow,
    WorkflowResult,
    molecule_to_dict,
    molecule_to_stjames,
    register_result,
)


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
        return f"<BDEResult energy={self.energy} H bdes={n}>"

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
    initial_molecule: MoleculeInput,
    mode: str = "rapid",
    atoms: list[int] | None = None,
    all_CH: bool = False,
    all_CX: bool = False,
    name: str = "BDE Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
    is_draft: bool = False,
) -> Workflow:
    """
    Submits a Bond-Dissociation Energy (BDE) workflow to the API.

    :param initial_molecule: Molecule to calculate BDEs for.
    :param mode: Mode to run the calculation in.
    :param atoms: Atom indices (1-indexed) to dissociate.
    :param all_CH: Whether to dissociate all C-H bonds.
    :param all_CX: Whether to dissociate all C-X bonds (X = halogen).
    :param name: Name of the workflow.
    :param folder_uuid: UUID of the folder to place the workflow in.
    :param folder: Folder object to store the workflow in.
    :param max_credits: Maximum number of credits to use for the workflow.
    :param webhook_url: URL that Rowan will POST to when the workflow completes.
    :param is_draft: If True, submit the workflow as a draft without starting execution.
    :returns: Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder:
        folder_uuid = folder.uuid
    initial_molecule = molecule_to_dict(initial_molecule)

    workflow = stjames.BDEWorkflow(
        initial_molecule=initial_molecule,
        mode=mode,
        atoms=atoms or [],
        all_CH=all_CH,
        all_CX=all_CX,
    )

    data = {
        "workflow_type": "bde",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": initial_molecule,
        "name": name,
        "folder_uuid": folder_uuid,
        "max_credits": max_credits,
        "webhook_url": webhook_url,
        "is_draft": is_draft,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


def find_ch_bonds(molecule: MoleculeInput, distance_max: float = 1.2) -> list[tuple[int, int]]:
    """
    Find all C-H bonds in a molecule.

    :param molecule: Molecule to search (Molecule, stjames.Molecule, or dict).
    :param distance_max: Maximum C-H distance to consider a bond (A).
    :returns: List of (carbon_index, hydrogen_index) tuples (1-based indices).

    Example::

        mol = Molecule.from_smiles("CCO")  # ethanol
        bonds = find_ch_bonds(mol)
        # [(1, 4), (1, 5), (1, 6), (2, 7), (2, 8)]
    """
    stj = molecule_to_stjames(molecule)
    return list(_find_CH_bonds(stj, distance_max))


def find_cx_bonds(molecule: MoleculeInput) -> list[tuple[int, int]]:
    """
    Find all C-X bonds in a molecule (X = F, Cl, Br, I, At, Ts).

    :param molecule: Molecule to search (Molecule, stjames.Molecule, or dict).
    :returns: List of (carbon_index, halogen_index) tuples (1-based indices).

    Example::

        mol = Molecule.from_smiles("CCCl")  # chloroethane
        bonds = find_cx_bonds(mol)
        # [(2, 3)]
    """
    stj = molecule_to_stjames(molecule)
    return list(_find_CX_bonds(stj))


def find_bonds(
    molecule: MoleculeInput,
    element_a: int,
    element_b: int,
    distance_max: float,
) -> list[tuple[int, int]]:
    """
    Find all bonds between two element types in a molecule.

    :param molecule: Molecule to search (Molecule, stjames.Molecule, or dict).
    :param element_a: Atomic number of first element.
    :param element_b: Atomic number of second element.
    :param distance_max: Maximum distance to consider a bond (A).
    :returns: List of (atom_a_index, atom_b_index) tuples (1-based indices).

    Example::

        mol = Molecule.from_smiles("O")  # water
        bonds = find_bonds(mol, 8, 1, 1.1)  # O-H bonds
        # [(1, 2), (1, 3)]
    """
    stj = molecule_to_stjames(molecule)
    return list(_find_AB_bonds(stj, element_a, element_b, distance_max))
