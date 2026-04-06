"""MacropKa workflow - predict macroscopic pKa values."""

from dataclasses import dataclass

import stjames

from ..folder import Folder
from ..types import MoleculeInput
from ..utils import api_client
from .base import Workflow, WorkflowResult, extract_smiles, register_result


@dataclass(frozen=True, slots=True)
class MacropKaMicrostate:
    """A microstate from a macropKa calculation."""

    smiles: str
    energy: float
    charge: int


@dataclass(frozen=True, slots=True)
class MacropKaValue:
    """A macroscopic pKa value."""

    initial_charge: int
    final_charge: int
    pka: float


@register_result("macropka")
class MacropKaResult(WorkflowResult):
    """Result from a macropKa workflow."""

    _stjames_class = stjames.MacropKaWorkflow

    def __repr__(self) -> str:
        iep = self.isoelectric_point
        n = len(self.pka_values)
        return f"<MacropKaResult isoelectric_point={iep} pka_values={n}>"

    @property
    def isoelectric_point(self) -> float | None:
        """Isoelectric point (pH units)."""
        return self._workflow.isoelectric_point

    @property
    def solvation_energy(self) -> float | None:
        """Solvation energy (kcal/mol)."""
        return self._workflow.solvation_energy

    @property
    def kpuu_probability(self) -> float | None:
        """Probability that Kpuu >= 0.3."""
        return self._workflow.kpuu_probability

    @property
    def microstates(self) -> list[MacropKaMicrostate]:
        """All microstates."""
        return [
            MacropKaMicrostate(
                smiles=m.smiles,
                energy=m.energy,
                charge=m.charge,
            )
            for m in self._workflow.microstates
        ]

    @property
    def pka_values(self) -> list[MacropKaValue]:
        """Macroscopic pKa values."""
        return [
            MacropKaValue(
                initial_charge=v.initial_charge,
                final_charge=v.final_charge,
                pka=v.pKa,
            )
            for v in self._workflow.pKa_values
        ]

    @property
    def microstate_weights_by_ph(self) -> list[tuple[float, list[float]]]:
        """Microstate weights by pH as (pH, weights) pairs.

        Each weights list corresponds to the microstates in the same order
        as the `microstates` property.
        """
        return list(self._workflow.microstate_weights_by_pH)

    @property
    def logd_by_ph(self) -> list[tuple[float, float]]:
        """Distribution constant by pH as (pH, logD) pairs."""
        return list(self._workflow.logD_by_pH)

    @property
    def aqueous_solubility_by_ph(self) -> list[tuple[float, float]]:
        """Aqueous solubility by pH as (pH, log(S)/L) pairs."""
        return list(self._workflow.aqueous_solubility_by_pH)


def submit_macropka_workflow(
    initial_smiles: str | MoleculeInput,
    min_pH: int = 0,
    max_pH: int = 14,
    min_charge: int = -2,
    max_charge: int = 2,
    compute_solvation_energy: bool = False,
    compute_aqueous_solubility: bool = True,
    name: str = "Macropka Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
    is_draft: bool = False,
) -> Workflow:
    """
    Submits a macropKa workflow to the API.

    :param initial_smiles: Molecule to calculate macroscopic pKa values for. Accepts
        a SMILES string or any molecule type (RowanMolecule, stjames.Molecule, RDKit Mol,
        or dict). The molecule must have a SMILES string associated with it, as macropKa
        models are 2D/SMILES-based and do not use 3D coordinates.
    :param min_pH: Minimum pH to use in the macropka workflow.
    :param max_pH: Maximum pH to use in the macropka workflow.
    :param min_charge: Minimum charge to use in the macropka workflow.
    :param max_charge: Maximum charge to use in the macropka workflow.
    :param compute_solvation_energy: Whether to compute the solvation energy.
    :param compute_aqueous_solubility: Whether to compute aqueous solubility for each pH.
    :param name: Name of the workflow.
    :param folder_uuid: UUID of the folder to store the workflow in.
    :param folder: Folder object to store the workflow in.
    :param max_credits: Maximum number of credits to use for the workflow.
    :param webhook_url: URL that Rowan will POST to when the workflow completes.
    :param is_draft: If True, submit the workflow as a draft without starting execution.
    :returns: Workflow object representing the submitted workflow.
    :raises ValueError: If the molecule has no SMILES associated with it.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder:
        folder_uuid = folder.uuid
    initial_smiles = extract_smiles(initial_smiles)
    workflow = stjames.MacropKaWorkflow(
        initial_smiles=initial_smiles,
        min_pH=min_pH,
        max_pH=max_pH,
        min_charge=min_charge,
        max_charge=max_charge,
        compute_solvation_energy=compute_solvation_energy,
        compute_aqueous_solubility=compute_aqueous_solubility,
    )

    data = {
        "workflow_type": "macropka",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_smiles": initial_smiles,
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
