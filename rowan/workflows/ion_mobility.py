"""Ion-mobility workflow - predict collision cross sections."""

import stjames

from ..folder import Folder
from ..types import MoleculeInput
from ..utils import api_client
from .base import Workflow, WorkflowResult, molecule_to_dict, register_result


@register_result("ion_mobility")
class IonMobilityResult(WorkflowResult):
    """Result from an ion-mobility workflow."""

    _stjames_class = stjames.IonMobilityWorkflow

    def __repr__(self) -> str:
        ccs = self.average_ccs
        std = self.average_ccs_stdev
        return f"<IonMobilityResult average_ccs={ccs} A^2 stdev={std} A^2>"

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
    initial_molecule: MoleculeInput,
    temperature: float = 300,
    protonate: bool = False,
    do_csearch: bool = True,
    do_optimization: bool = True,
    name: str = "Ion-Mobility Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
) -> Workflow:
    """
    Submits an ion-mobility workflow to the API.

    :param initial_molecule: Molecule used in the scan.
    :param temperature: Temperature at which to predict CCS values (K).
    :param protonate: Whether or not to automatically detect protonation site.
        If `True`, every basic site will be protonated and values returned for the most stable.
    :param do_csearch: Whether to perform a conformational search on the molecule.
    :param do_optimization: Whether to perform an optimization on the molecule.
    :param name: Name of the workflow.
    :param folder_uuid: UUID of the folder to store the workflow in.
    :param folder: Folder object to store the workflow in.
    :param max_credits: Maximum number of credits to use for the workflow.
    :param webhook_url: URL that Rowan will POST to when the workflow completes.
    :returns: Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder:
        folder_uuid = folder.uuid
    initial_molecule = molecule_to_dict(initial_molecule)

    workflow = stjames.IonMobilityWorkflow(
        initial_molecule=initial_molecule,
        temperature=temperature,
        protonate=protonate,
        do_csearch=do_csearch,
        do_optimization=do_optimization,
    )

    data = {
        "workflow_type": "ion_mobility",
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
