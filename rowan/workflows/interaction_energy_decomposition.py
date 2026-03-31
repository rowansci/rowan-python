"""Interaction energy decomposition workflow - SAPT0 decomposition between molecular fragments."""

from typing import Literal

import stjames
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

from ..folder import Folder
from ..types import MoleculeInput
from ..utils import api_client
from .base import Workflow, WorkflowResult, molecule_to_dict, molecule_to_stjames, register_result


def _validate_fragment_separation(molecule: MoleculeInput, fragment1_indices: list[int]) -> None:
    """
    Validate that no atom in fragment 1 is covalently bonded to an atom outside fragment 1.

    :raises ValueError: If a cross-fragment covalent bond is detected.
    """
    stj_mol = molecule_to_stjames(molecule)
    xyz = stj_mol.to_xyz()

    rdmol = Chem.MolFromXYZBlock(xyz)  # type: ignore[attr-defined]
    if rdmol is None:
        return  # Can't determine connectivity, skip validation
    rdDetermineBonds.DetermineConnectivity(rdmol)

    frag1 = {i - 1 for i in fragment1_indices}  # convert to 0-indexed

    for bond in rdmol.GetBonds():
        a = bond.GetBeginAtomIdx()
        b = bond.GetEndAtomIdx()
        if (a in frag1) != (b in frag1):
            raise ValueError(
                f"Atom {a + 1} and atom {b + 1} are covalently bonded but assigned to different "
                f"fragments. SAPT requires non-covalently bound fragments."
            )


@register_result("interaction_energy_decomposition")
class InteractionEnergyDecompositionResult(WorkflowResult):
    """Result from an interaction energy decomposition workflow."""

    _stjames_class = stjames.InteractionEnergyDecompositionWorkflow

    def __repr__(self) -> str:
        e = self.total_interaction_energy
        return f"<InteractionEnergyDecompositionResult total_interaction_energy={e}>"

    @property
    def fragment1_indices(self) -> list[int]:
        """Atom indices (1-indexed) defining fragment 1."""
        return list(self._workflow.fragment1_indices)

    @property
    def total_interaction_energy(self) -> float | None:
        """Total interaction energy (kcal/mol)."""
        r = self._workflow.energy_decomposition_result
        return r.total_interaction_energy if r else None

    @property
    def electrostatic_interaction_energy(self) -> float | None:
        """Electrostatic interaction energy (kcal/mol)."""
        r = self._workflow.energy_decomposition_result
        return r.electrostatic_interaction_energy if r else None

    @property
    def exchange_interaction_energy(self) -> float | None:
        """Exchange interaction energy (kcal/mol)."""
        r = self._workflow.energy_decomposition_result
        return r.exchange_interaction_energy if r else None

    @property
    def dispersion_interaction_energy(self) -> float | None:
        """Dispersion interaction energy (kcal/mol)."""
        r = self._workflow.energy_decomposition_result
        return r.dispersion_interaction_energy if r else None

    @property
    def induction_interaction_energy(self) -> float | None:
        """Induction interaction energy (kcal/mol)."""
        r = self._workflow.energy_decomposition_result
        return r.induction_interaction_energy if r else None


def submit_interaction_energy_decomposition_workflow(
    initial_molecule: MoleculeInput,
    fragment1_indices: list[int],
    method: Literal["sapt0"] = "sapt0",
    basis_set: str = "jun-cc-pVDZ",
    name: str = "Interaction Energy Decomposition",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
) -> Workflow:
    """
    Submits an interaction energy decomposition (SAPT0) workflow to the API.

    Decomposes the interaction energy between two molecular fragments into
    electrostatic, exchange, dispersion, and induction components.

    :param initial_molecule: Dimer molecule to decompose.
    :param fragment1_indices: Atom indices (1-indexed) defining fragment 1.
        Fragment 2 is all remaining atoms.
    :param method: Energy decomposition method. Currently only ``"sapt0"`` is supported.
    :param basis_set: Basis set for the calculation. Defaults to ``"jun-cc-pVDZ"``.
    :param name: Name of the workflow.
    :param folder_uuid: UUID of the folder to place the workflow in.
    :param folder: Folder object to store the workflow in.
    :param max_credits: Maximum number of credits to use for the workflow.
    :param webhook_url: URL that Rowan will POST to when the workflow completes.
    :returns: Workflow object representing the submitted workflow.
    :raises ValueError: If both folder and folder_uuid are provided, or if fragment 1
        contains atoms covalently bonded to atoms outside fragment 1.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder:
        folder_uuid = folder.uuid

    _validate_fragment_separation(initial_molecule, fragment1_indices)

    initial_molecule = molecule_to_dict(initial_molecule)

    energy_decomposition_settings = stjames.EnergyDecompositionSettings(
        method=method,
        basis_set=stjames.BasisSet(name=basis_set),
    )

    workflow = stjames.InteractionEnergyDecompositionWorkflow(
        initial_molecule=initial_molecule,
        fragment1_indices=fragment1_indices,
        energy_decomposition_settings=energy_decomposition_settings,
    )

    data = {
        "workflow_type": "interaction_energy_decomposition",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": initial_molecule,
        "name": name,
        "folder_uuid": folder_uuid,
        "max_credits": max_credits,
        "webhook_url": webhook_url,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())
