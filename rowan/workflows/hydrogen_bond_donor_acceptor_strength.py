"""Hydrogen-bond donor/acceptor-strength workflow."""

import math
from dataclasses import dataclass

import stjames

from ..folder import Folder
from ..types import MoleculeInput
from ..utils import api_client
from .base import Workflow, WorkflowResult, molecule_to_dict, register_result


@dataclass(frozen=True, slots=True)
class HydrogenBondAcceptorSite:
    """Hydrogen-bond-acceptor site."""

    atom_idx: int
    pkbhx: float
    position: tuple[float, float, float]
    name: str | None = None


@dataclass(frozen=True, slots=True)
class HydrogenBondDonorSite:
    """Hydrogen-bond-donor site."""

    atom_idx: int
    pk_alpha: float
    position: tuple[float, float, float]


@register_result("hydrogen_bond_basicity")
class HydrogenBondDonorAcceptorStrengthResult(WorkflowResult):
    """Result from a hydrogen-bond donor/acceptor-strength workflow."""

    _stjames_class = stjames.HydrogenBondBasicityWorkflow

    def __repr__(self) -> str:
        n_acc = len(self.acceptor_sites)
        n_don = len(self.donor_sites)
        return f"<HydrogenBondDonorAcceptorStrengthResult acceptors={n_acc} donors={n_don}>"

    @property
    def acceptor_sites(self) -> list[HydrogenBondAcceptorSite]:
        """Hydrogen bond acceptor sites with pKBHX values."""
        return [
            HydrogenBondAcceptorSite(
                atom_idx=s.atom_idx,
                pkbhx=s.pkbhx,
                position=s.position,
                name=s.name,
            )
            for s in self._workflow.hba_sites
        ]

    @property
    def donor_sites(self) -> list[HydrogenBondDonorSite]:
        """Hydrogen bond donor sites with pK_alpha values."""
        return [
            HydrogenBondDonorSite(
                atom_idx=s.atom_idx,
                pk_alpha=s.pk_alpha,
                position=s.position,
            )
            for s in self._workflow.hbd_sites
        ]

    @property
    def molecular_pkbhx(self) -> float | None:
        """Overall molecular HBA strength as log10(sum of 10^pkbhx) for sites with pkbhx > -1."""
        sites = [s.pkbhx for s in self.acceptor_sites if s.pkbhx > -1]
        if not sites:
            return None
        return math.log10(sum(10**pkbhx for pkbhx in sites))

    @property
    def molecular_pk_alpha(self) -> float | None:
        """Overall molecular HBD strength as log10(sum of 10^pk_alpha) (pk_alpha > -1)."""
        sites = [s.pk_alpha for s in self.donor_sites if s.pk_alpha > -1]
        if not sites:
            return None
        return math.log10(sum(10**pk_alpha for pk_alpha in sites))


def submit_hydrogen_bond_donor_acceptor_strength_workflow(
    initial_molecule: MoleculeInput,
    do_csearch: bool = True,
    do_optimization: bool = True,
    name: str = "Hydrogen-Bond Acceptor/Donor Strength",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
    is_draft: bool = False,
) -> Workflow:
    """
    Submits a hydrogen-bond donor/acceptor-strength workflow to the API.

    :param initial_molecule: Molecule to calculate HBA/HBD strength for.
    :param do_csearch: Whether to perform a conformational search.
    :param do_optimization: Whether to perform an optimization.
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

    workflow = stjames.HydrogenBondBasicityWorkflow(
        initial_molecule=initial_molecule,
        do_csearch=do_csearch,
        do_optimization=do_optimization,
    )

    data = {
        "workflow_type": "hydrogen_bond_basicity",
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


# Backwards compatibility aliases
HydrogenBondBasicityResult = HydrogenBondDonorAcceptorStrengthResult
submit_hydrogen_bond_basicity_workflow = submit_hydrogen_bond_donor_acceptor_strength_workflow
