"""Ion mobility workflow - predict collision cross sections."""

from typing import Any

import stjames
from rdkit import Chem

from ..utils import api_client
from .base import RdkitMol, StJamesMolecule, Workflow, WorkflowResult, register_result


@register_result("ion_mobility")
class IonMobilityResult(WorkflowResult):
    """Result from an ion mobility workflow."""

    _stjames_class = stjames.IonMobilityWorkflow

    def __repr__(self) -> str:
        ccs = self.average_ccs
        std = self.average_ccs_stdev
        return f"<IonMobilityResult average_ccs={ccs} stdev={std}>"

    @property
    def average_ccs(self) -> float | None:
        """Average collision cross section (Angstrom^2)."""
        return self._workflow.average_ccs

    @property
    def average_ccs_stdev(self) -> float | None:
        """Uncertainty in average CCS."""
        return self._workflow.average_ccs_stdev

    @property
    def conformer_ccs(self) -> list[float]:
        """Collision cross section per conformer (Angstrom^2)."""
        return list(self._workflow.conformer_ccs)

    @property
    def boltzmann_weights(self) -> list[float]:
        """Boltzmann weights for conformers."""
        return list(self._workflow.boltzmann_weights)


def submit_ion_mobility_workflow(
    initial_molecule: dict[str, Any] | StJamesMolecule | RdkitMol,
    temperature: float = 300,
    protonate: bool = False,
    do_csearch: bool = True,
    do_optimization: bool = True,
    name: str = "Ion-Mobility Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits an ion-mobility workflow to the API.

    :param initial_molecule: The molecule used in the scan.
    :param temperature: The temperature at which to predict CCS values.
    :param protonate: Whether or not to automatically detect protonation site.
        If `True`, every basic site will be protonated and values returned for the most stable.
    :param do_csearch: Whether to perform a conformational search on the molecule.
    :param do_optimization: Whether to perform an optimization on the molecule.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to store the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if isinstance(initial_molecule, StJamesMolecule):
        initial_molecule = initial_molecule.model_dump(mode="json")
    elif isinstance(initial_molecule, Chem.rdchem.Mol | Chem.rdchem.RWMol):
        initial_molecule = StJamesMolecule.from_rdkit(initial_molecule, cid=0)

    workflow = stjames.IonMobilityWorkflow(
        initial_molecule=initial_molecule,
        temperature=temperature,
        protonate=protonate,
        do_csearch=do_csearch,
        do_optimization=do_optimization,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "ion_mobility",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": initial_molecule,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


__all__ = ["IonMobilityResult", "submit_ion_mobility_workflow"]
