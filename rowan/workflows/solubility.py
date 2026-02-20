"""Solubility workflow - predict molecular solubility in various solvents."""

from dataclasses import dataclass
from typing import Literal

import stjames

from ..utils import api_client
from .base import Workflow, WorkflowResult, register_result


@dataclass(frozen=True, slots=True)
class SolubilityEntry:
    """Solubility results for a single solvent."""

    solvent: str
    solubilities: tuple[float, ...]
    uncertainties: tuple[float | None, ...]


@register_result("solubility")
class SolubilityResult(WorkflowResult):
    """Result from an aqueous solubility workflow."""

    _stjames_class = stjames.SolubilityWorkflow

    def __repr__(self) -> str:
        solvents = [s.solvent for s in self.solubilities]
        return f"<SolubilityResult solvents={solvents}>"

    @property
    def solubilities(self) -> list[SolubilityEntry]:
        """Solubility results per solvent."""
        return [
            SolubilityEntry(
                solvent=solvent,
                solubilities=tuple(result.solubilities),
                uncertainties=tuple(result.uncertainties),
            )
            for solvent, result in self._workflow.solubilities.items()
        ]

    @property
    def temperatures(self) -> list[float]:
        """Temperatures in Kelvin."""
        return list(self._workflow.temperatures)


def submit_solubility_workflow(
    initial_smiles: str,
    solubility_method: Literal["fastsolv", "kingfisher", "esol"] = "fastsolv",
    solvents: list[str] | None = None,
    temperatures: list[float] | None = None,
    name: str = "Solubility Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a solubility workflow to the API.

    :param initial_smiles: The smiles of the molecule to calculate the solubility of.
    :param solubility_method: The name of the desired model for solubility prediction.
    :param solvents: The list of solvents to use for the calculation.
    :param temperatures: The list of temperatures to use for the calculation (K).
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if not solvents:
        solvents = ["CCCCCC", "CC1=CC=CC=C1", "C1CCCO1", "CC(=O)OCC", "CCO", "CC#N"]

    if not temperatures:
        temperatures = [273.15, 298.15, 323.15, 348.15, 373.15]

    workflow = stjames.SolubilityWorkflow(
        initial_smiles=initial_smiles,
        solubility_method=solubility_method,
        solvents=solvents,
        temperatures=temperatures,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "solubility",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_smiles": initial_smiles,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


__all__ = ["SolubilityEntry", "SolubilityResult", "submit_solubility_workflow"]
