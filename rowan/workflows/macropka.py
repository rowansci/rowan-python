"""MacropKa workflow - predict macroscopic pKa values."""

from dataclasses import dataclass

import stjames

from ..utils import api_client
from .base import Workflow, WorkflowResult, register_result


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
    def logd_by_ph(self) -> list[tuple[float, float]]:
        """Distribution constant by pH as (pH, logD) pairs."""
        return list(self._workflow.logD_by_pH)

    @property
    def aqueous_solubility_by_ph(self) -> list[tuple[float, float]]:
        """Aqueous solubility by pH as (pH, log(S)/L) pairs."""
        return list(self._workflow.aqueous_solubility_by_pH)


def submit_macropka_workflow(
    initial_smiles: str,
    min_pH: int = 0,
    max_pH: int = 14,
    min_charge: int = -2,
    max_charge: int = 2,
    compute_solvation_energy: bool = False,
    compute_aqueous_solubility: bool = False,
    name: str = "Macropka Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a macropka workflow to the API.

    :param initial_smiles: The molecule used in the macropka workflow.
    :param min_pH: The minimum pH to use in the macropka workflow.
    :param max_pH: The maximum pH to use in the macropka workflow.
    :param min_charge: The minimum charge to use in the macropka workflow.
    :param max_charge: The maximum charge to use in the macropka workflow.
    :param compute_solvation_energy: Whether to compute the solvation energy.
    :param compute_aqueous_solubility: Whether to compute the aqueous solubility for each pH.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to store the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
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
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "macropka",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_smiles": initial_smiles,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


__all__ = [
    "MacropKaMicrostate",
    "MacropKaResult",
    "MacropKaValue",
    "submit_macropka_workflow",
]
