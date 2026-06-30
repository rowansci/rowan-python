"""Binding affinity workflow — SQM-based scoring of protein–ligand complexes."""

from dataclasses import dataclass

import stjames
from stjames import SinglePointEnergySettings

from ..folder import Folder
from ..protein import Protein
from ..types import StructureInput
from ..utils import api_client
from .base import (
    Message,
    Workflow,
    WorkflowResult,
    molecule_to_dict,
    parse_messages,
    register_result,
    require_coordinates,
)


@dataclass(frozen=True, slots=True)
class BindingAffinityScore:
    """
    Binding affinity score for a single pose.

    :param binding_affinity: binding affinity in kcal/mol
    :param strain: strain energy in kcal/mol, or None if not computed
    """

    binding_affinity: float
    strain: float | None


@register_result("binding_affinity")
class BindingAffinityResult(WorkflowResult):
    """Result from a binding affinity workflow."""

    _stjames_class = stjames.BindingAffinityWorkflow

    def __repr__(self) -> str:
        n = len(self.scores)
        return f"<BindingAffinityResult scores={n}>"

    @property
    def scores(self) -> list[BindingAffinityScore]:
        """Binding affinity scores for each scored pose."""
        return [
            BindingAffinityScore(
                binding_affinity=r.binding_affinity,
                strain=r.strain,
            )
            for r in (self._workflow.binding_affinity_results or [])
        ]

    @property
    def messages(self) -> list[Message]:
        """Any messages or warnings from the workflow."""
        return parse_messages(self._workflow.messages)


def submit_binding_affinity_workflow(
    protein: str | Protein,
    ligand_residue_name: str | None = None,
    ligand_structures: list[StructureInput] | None = None,
    binding_affinity_settings: SinglePointEnergySettings | None = None,
    name: str = "Binding Affinity Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
    is_draft: bool = False,
) -> Workflow:
    """
    Submits a binding affinity workflow to the API.

    Scores ligand poses using SQM-based energies. Two submission modes:

    **Mode 1 — holo protein:** protein already contains the bound ligand. Pass
    ``ligand_residue_name`` to identify which residue is the ligand vs. the receptor.
    Do not pass ``ligand_structures``.

    **Mode 2 — apo protein + external poses:** protein has no bound ligand. Pass
    ``ligand_structures`` with poses that are already in the protein's coordinate frame.
    Do not pass ``ligand_residue_name``. When scoring multiple ligands, prefer this mode
    over separate per-ligand workflows — all poses share the same pocket geometry.

    :param protein: protein structure. Can be input as a UUID or a Protein object.
    :param ligand_residue_name: residue name identifying the ligand in a holo protein PDB
        (mode 1 only).
    :param ligand_structures: external ligand poses to score, already in the protein's
        coordinate frame. Must have 3D coordinates (mode 2 only).
    :param binding_affinity_settings: SQM settings controlling geometry optimization and
        energy evaluation. Defaults to PM6-D3H4X/COSMO optimization followed by
        PM6-D3H4X/COSMO2 single-point in water.
    :param name: name of the workflow.
    :param folder_uuid: UUID of the folder to place the workflow in.
    :param folder: Folder object to store the workflow in.
    :param max_credits: maximum number of credits to use for the workflow.
    :param webhook_url: URL that Rowan will POST to when the workflow completes.
    :param is_draft: if True, submit the workflow as a draft without starting execution.
    :returns: Workflow object representing the submitted workflow.
    :raises ValueError: if folder arguments conflict.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder:
        folder_uuid = folder.uuid
    if isinstance(protein, Protein):
        protein = protein.uuid

    mol_dicts: list[dict] = []
    if ligand_structures:
        for mol in ligand_structures:
            require_coordinates(mol)
            mol_dicts.append(molecule_to_dict(mol))

    workflow = stjames.BindingAffinityWorkflow(
        protein=protein,
        ligand_residue_name=ligand_residue_name,
        ligand_structures=mol_dicts or [],
        binding_affinity_settings=binding_affinity_settings or SinglePointEnergySettings(),
    )

    data = {
        "workflow_type": "binding_affinity",
        "workflow_data": workflow.model_dump(serialize_as_any=True, mode="json"),
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
