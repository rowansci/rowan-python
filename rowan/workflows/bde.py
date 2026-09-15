"""BDE workflow - Bond Dissociation Energy calculations."""

from dataclasses import dataclass
from typing import Any

import stjames
from stjames import MultiStageOptSettings
from stjames.workflows.bde import find_AB_bonds as _find_AB_bonds
from stjames.workflows.bde import find_CH_bonds as _find_CH_bonds
from stjames.workflows.bde import find_CX_bonds as _find_CX_bonds

from ..folder import Folder
from ..utils import api_client
from .base import (
    StructureInput,
    Workflow,
    WorkflowResult,
    molecule_to_dict,
    molecule_to_stjames,
    register_result,
    require_coordinates,
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
    initial_molecule: StructureInput,
    mode: str = "omol25_conserving_s",
    multistage_opt_settings: MultiStageOptSettings | None = None,
    fragment_indices: list[list[int]] | None = None,
    all_CH: bool = False,
    all_CX: bool = False,
    name: str = "BDE Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
    is_draft: bool = False,
) -> Workflow:
    """Submits a Bond-Dissociation Energy (BDE) workflow to the API.

    Args:
        initial_molecule: molecule to calculate BDEs for
        mode: level of theory to run the calculation at, given as a method string:
            - `omol25_conserving_s` – neural network potential (default)
            - `g_xtb//gfn2_xtb` – semiempirical
            - `r2scan3c//gfn2_xtb` – DFT single point on a semiempirical geometry
        multistage_opt_settings: explicit method sequence to use instead of the one `mode` would
            pick – the optimization stage(s) followed by a final singlepoint, given as a
            `MultiStageOptSettings`. When omitted, the sequence is built automatically from `mode`.
            When
            supplied, it replaces that sequence
        fragment_indices: 1-indexed atoms of each fragment to dissociate. Each fragment must
            connect to the rest of the molecule by a single bond
        all_CH: whether to dissociate all C-H bonds
        all_CX: whether to dissociate all C-X bonds (X = halogen)
        name: name of the workflow
        folder_uuid: UUID of the folder to place the workflow in
        folder: destination folder
        max_credits: maximum credits for the workflow
        webhook_url: URL that Rowan will POST to when the workflow completes
        is_draft: save as a draft without starting execution

    Returns:
        submitted workflow

    Raises:
        httpx.HTTPStatusError: request to the API fails
    """
    require_coordinates(initial_molecule)
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder:
        folder_uuid = folder.uuid
    mol_dict = molecule_to_dict(initial_molecule)

    workflow_kwargs: dict[str, Any] = {
        "initial_molecule": mol_dict,
        "mode": mode,
        "fragment_indices": fragment_indices or [],
        "all_CH": all_CH,
        "all_CX": all_CX,
    }
    if multistage_opt_settings is not None:
        workflow_kwargs["multistage_opt_settings"] = multistage_opt_settings

    workflow = stjames.BDEWorkflow.model_validate(workflow_kwargs)

    data = {
        "workflow_type": "bde",
        "workflow_data": workflow.model_dump(mode="json", serialize_as_any=True),
        "initial_molecule": mol_dict,
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


def find_ch_bonds(molecule: StructureInput, distance_max: float = 1.2) -> list[tuple[int, int]]:
    """Find all C-H bonds in a molecule.

    Args:
        molecule: molecule to search (Molecule, stjames.Molecule, or dict)
        distance_max: maximum C-H distance to consider a bond (A)

    Returns:
        list of (carbon_index, hydrogen_index) tuples (1-based indices)

    Examples:
        ```python
        mol = Molecule.from_smiles("CCO")  # ethanol
        bonds = find_ch_bonds(mol)
        # [(1, 4), (1, 5), (1, 6), (2, 7), (2, 8)]
        ```
    """
    stj = molecule_to_stjames(molecule)
    return list(_find_CH_bonds(stj, distance_max))


def find_cx_bonds(molecule: StructureInput) -> list[tuple[int, int]]:
    """Find all C-X bonds in a molecule (X = F, Cl, Br, I, At, Ts).

    Args:
        molecule: molecule to search (Molecule, stjames.Molecule, or dict)

    Returns:
        list of (carbon_index, halogen_index) tuples (1-based indices)

    Examples:
        ```python
        mol = Molecule.from_smiles("CCCl")  # chloroethane
        bonds = find_cx_bonds(mol)
        # [(2, 3)]
        ```
    """
    stj = molecule_to_stjames(molecule)
    return list(_find_CX_bonds(stj))


def find_bonds(
    molecule: StructureInput,
    element_a: int,
    element_b: int,
    distance_max: float,
) -> list[tuple[int, int]]:
    """Find all bonds between two element types in a molecule.

    Args:
        molecule: molecule to search (Molecule, stjames.Molecule, or dict)
        element_a: atomic number of first element
        element_b: atomic number of second element
        distance_max: maximum distance to consider a bond (A)

    Returns:
        list of (atom_a_index, atom_b_index) tuples (1-based indices)

    Examples:
        ```python
        mol = Molecule.from_smiles("O")  # water
        bonds = find_bonds(mol, 8, 1, 1.1)  # O-H bonds
        # [(1, 2), (1, 3)]
        ```

    Same-element searches return unique undirected bonds without self-pairs:

        >>> peroxide = stjames.Molecule.from_smiles("OO")
        >>> find_bonds(peroxide, 8, 8, 1.7)
        [(1, 2)]
    """
    stj = molecule_to_stjames(molecule)
    bonds = _find_AB_bonds(stj, element_a, element_b, distance_max)
    if element_a != element_b:
        return list(bonds)
    return sorted(
        {(min(atom_a, atom_b), max(atom_a, atom_b)) for atom_a, atom_b in bonds if atom_a != atom_b}
    )
