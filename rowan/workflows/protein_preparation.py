"""Protein preparation workflow - prepare protein structures for simulation."""

from typing import Literal

import stjames

from ..folder import Folder
from ..protein import Protein, retrieve_protein
from ..types import ProteinUUID
from ..utils import api_client
from .base import Workflow, WorkflowResult, register_result

_DEFAULT_RETAIN_NON_POLYMER: dict[str | int, str | None] = {
    "NA": None,
    "CL": None,
    "MG": None,
}


@register_result("protein_preparation")
class ProteinPreparationResult(WorkflowResult):
    """Result from a protein-preparation workflow."""

    _stjames_class = stjames.ProteinPreparationWorkflow

    def __repr__(self) -> str:
        return f"<ProteinPreparationResult prepared_protein_uuid={self.prepared_protein_uuid}>"

    @property
    def prepared_protein_uuid(self) -> ProteinUUID | None:
        """UUID of the prepared protein structure."""
        return getattr(self._workflow, "prepared_protein", None)

    def get_prepared_protein(self) -> Protein:
        """
        Fetch the prepared protein structure.

        .. note::
            Makes one API call on first access.
            Results are cached. Call clear_cache() to refresh.

        :returns: prepared Protein object
        :raises ValueError: if the workflow has not produced a prepared protein
        """
        if not (uuid := self.prepared_protein_uuid):
            raise ValueError("Protein preparation has no prepared protein UUID")
        if "prepared_protein" not in self._cache:
            protein = retrieve_protein(uuid, workflow_uuid=self.workflow_uuid)
            if protein.data is None:
                protein.refresh()
            if protein.data is None:
                raise ValueError("Prepared protein record has no structure data")
            self._cache["prepared_protein"] = protein
        return self._cache["prepared_protein"]


def submit_protein_preparation_workflow(
    protein: Protein | ProteinUUID,
    add_missing_method: Literal["boltz_2", "pdbfixer"] | None = "boltz_2",
    cap_residues: Literal["ace_nme", "terminal_templates"] | None = "ace_nme",
    protonation_method: Literal["openmm", "protonate_utils", "propka_3"] = "openmm",
    pH: float = 7.4,
    retain_protonation: bool = True,
    retain_non_polymer: dict[str | int, str | None] | None = _DEFAULT_RETAIN_NON_POLYMER,
    name: str = "Protein Preparation Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
    is_draft: bool = False,
) -> Workflow:
    """
    Submit a protein-preparation workflow to the API.

    Full protein preparation can take around ten minutes, depending on the structure and
    settings. For a faster PDBFixer/OpenMM-only path, use ``Protein.prepare()``.

    :param protein: protein to prepare, as a UUID or Protein object
    :param add_missing_method: method for adding missing atoms and residues before protonation;
        None skips this step
    :param cap_residues: terminal-residue capping method; ``"ace_nme"`` requires an
        add-missing method and is incompatible with ``"protonate_utils"``; None disables capping
    :param protonation_method: method for adding hydrogens
    :param pH: pH used to determine protonation states
    :param retain_protonation: whether to retain existing protonation states
    :param retain_non_polymer: non-polymer residues to retain, keyed by residue name or 0-based
        residue index. Values are SMILES strings used for parameterization. Known ions and waters
        may map to None; other residues require a SMILES string. None removes all non-polymer
        residues.
    :param name: name of the workflow
    :param folder_uuid: UUID of the folder to place the workflow in
    :param folder: Folder object to store the workflow in
    :param max_credits: maximum number of credits to use for the workflow
    :param webhook_url: URL that Rowan will POST to when the workflow completes
    :param is_draft: if True, submit the workflow as a draft without starting execution
    :returns: Workflow object representing the submitted workflow
    :raises ValueError: if settings are incompatible or a retained residue mapping is invalid
    :raises requests.HTTPError: if the request to the API fails
    """
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder:
        folder_uuid = folder.uuid
    if isinstance(protein, Protein):
        protein = protein.uuid

    workflow = stjames.ProteinPreparationWorkflow(
        protein=protein,
        add_missing_method=add_missing_method,
        cap_residues=cap_residues,
        protonation_method=protonation_method,
        pH=pH,
        retain_protonation=retain_protonation,
        retain_non_polymer=(dict(retain_non_polymer) if retain_non_polymer is not None else None),
    )

    data = {
        "workflow_type": "protein_preparation",
        "workflow_data": workflow.model_dump(mode="json"),
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
