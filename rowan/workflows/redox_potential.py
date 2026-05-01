"""Redox potential workflow - calculate oxidation/reduction potentials."""

import stjames

from ..calculation import Calculation, retrieve_calculation
from ..folder import Folder
from ..utils import api_client
from .base import (
    Message,
    Mode,
    MoleculeInput,
    Workflow,
    WorkflowResult,
    molecule_to_dict,
    parse_messages,
    register_result,
)


@register_result("redox_potential")
class RedoxPotentialResult(WorkflowResult):
    """Result from a redox-potential workflow."""

    _stjames_class = stjames.RedoxPotentialWorkflow

    def __repr__(self) -> str:
        ox = self.oxidation_potential
        red = self.reduction_potential
        parts = []
        if ox is not None:
            parts.append(f"oxidation={ox:.3f}V")
        if red is not None:
            parts.append(f"reduction={red:.3f}V")
        return f"<RedoxPotentialResult {' '.join(parts)}>"

    @property
    def oxidation_potential(self) -> float | None:
        """Oxidation potential in V (vs SHE)."""
        return getattr(self._workflow, "oxidation_potential", None)

    @property
    def reduction_potential(self) -> float | None:
        """Reduction potential in V (vs SHE)."""
        return getattr(self._workflow, "reduction_potential", None)

    @property
    def neutral_molecule_uuid(self) -> str | None:
        """UUID of the optimized neutral molecule calculation."""
        return getattr(self._workflow, "neutral_molecule", None)

    @property
    def cation_molecule_uuid(self) -> str | None:
        """UUID of the optimized cation (oxidized) molecule calculation."""
        return getattr(self._workflow, "cation_molecule", None)

    @property
    def anion_molecule_uuid(self) -> str | None:
        """UUID of the optimized anion (reduced) molecule calculation."""
        return getattr(self._workflow, "anion_molecule", None)

    def get_neutral_molecule(self) -> Calculation | None:
        """Fetch the optimized neutral molecule calculation."""
        if not (uuid := self.neutral_molecule_uuid):
            return None
        if "neutral_molecule" not in self._cache:
            self._cache["neutral_molecule"] = retrieve_calculation(uuid)
        return self._cache["neutral_molecule"]

    def get_cation_molecule(self) -> Calculation | None:
        """Fetch the optimized cation (oxidized) molecule calculation."""
        if not (uuid := self.cation_molecule_uuid):
            return None
        if "cation_molecule" not in self._cache:
            self._cache["cation_molecule"] = retrieve_calculation(uuid)
        return self._cache["cation_molecule"]

    def get_anion_molecule(self) -> Calculation | None:
        """Fetch the optimized anion (reduced) molecule calculation."""
        if not (uuid := self.anion_molecule_uuid):
            return None
        if "anion_molecule" not in self._cache:
            self._cache["anion_molecule"] = retrieve_calculation(uuid)
        return self._cache["anion_molecule"]

    @property
    def messages(self) -> list[Message]:
        """Any messages or warnings from the workflow."""
        return parse_messages(getattr(self._workflow, "messages", None))


def submit_redox_potential_workflow(
    initial_molecule: MoleculeInput,
    reduction: bool = False,
    oxidation: bool = True,
    mode: Mode = Mode.RAPID,
    name: str = "Redox Potential Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
    is_draft: bool = False,
) -> Workflow:
    """
    Submits a redox-potential workflow to the API.

    :param initial_molecule: Molecule to calculate the redox potential of.
    :param reduction: Whether to calculate the reduction potential.
    :param oxidation: Whether to calculate the oxidation potential.
    :param mode: Mode to run the calculation in.
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

    workflow = stjames.RedoxPotentialWorkflow(
        initial_molecule=initial_molecule,
        oxidation=oxidation,
        reduction=reduction,
        mode=mode,
    )

    data = {
        "workflow_type": "redox_potential",
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
