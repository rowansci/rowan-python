"""Ion-mobility workflow - predict collision cross sections."""

import stjames

from rowan.folder import Folder
from rowan.utils import api_client

from .base import (
    StructureInput,
    Workflow,
    WorkflowResult,
    molecule_to_dict,
    register_result,
    require_coordinates,
)


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
    initial_molecule: StructureInput,
    temperature: float = 300,
    protonate: bool = False,
    do_csearch: bool = True,
    do_optimization: bool = True,
    name: str = "Ion-Mobility Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
    is_draft: bool = False,
) -> Workflow[IonMobilityResult]:
    """Submits an ion-mobility workflow to the API.

    Args:
        initial_molecule: molecule used in the scan
        temperature: temperature at which to predict CCS values (K)
        protonate: protonate each basic site and return values for the most stable form
        do_csearch: whether to perform a conformational search on the molecule.
            Requires do_optimization
        do_optimization: whether to perform an optimization on the molecule
        name: name of the workflow
        folder_uuid: UUID of the folder to store the workflow in
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
    if do_csearch and not do_optimization:
        raise ValueError(
            "`do_optimization` must be True when `do_csearch` is True; the conformers from "
            "the search must be optimized before the CCS calculation."
        )
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder:
        folder_uuid = folder.uuid
    mol_dict = molecule_to_dict(initial_molecule)

    workflow = stjames.IonMobilityWorkflow(
        initial_molecule=mol_dict,
        temperature=temperature,
        protonate=protonate,
        do_csearch=do_csearch,
        do_optimization=do_optimization,
    )

    data = {
        "workflow_type": "ion_mobility",
        "workflow_data": workflow.model_dump(mode="json"),
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
        return Workflow[IonMobilityResult](**response.json())
