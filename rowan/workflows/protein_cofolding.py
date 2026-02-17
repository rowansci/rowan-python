"""Protein cofolding workflow - predict protein-protein and protein-ligand complexes."""

from dataclasses import dataclass

import stjames

from ..utils import api_client
from .base import Workflow, WorkflowResult, register_result


@dataclass(frozen=True, slots=True)
class CofoldingResult:
    """A protein cofolding result."""

    average_plddt: float | None = None
    ptm: float | None = None
    iptm: float | None = None


@register_result("protein_cofolding")
class ProteinCofoldingResult(WorkflowResult):
    """Result from a protein cofolding workflow."""

    _stjames_class = stjames.ProteinCofoldingWorkflow

    def __repr__(self) -> str:
        results = self.cofolding_results
        n = len(results)
        if results:
            best = max(results, key=lambda r: r.iptm or 0)
            return f"<ProteinCofoldingResult results={n} best_iptm={best.iptm}>"
        return f"<ProteinCofoldingResult results={n}>"

    @property
    def cofolding_results(self) -> list[CofoldingResult]:
        """Cofolding results."""
        return [
            CofoldingResult(
                average_plddt=getattr(r, "avg_lddt", None),
                ptm=getattr(r, "ptm", None),
                iptm=getattr(r, "iptm", None),
            )
            for r in self._workflow.cofolding_results
        ]


def submit_protein_cofolding_workflow(
    initial_protein_sequences: list[str],
    initial_smiles_list: list[str] | None = None,
    ligand_binding_affinity_index: int | None = None,
    use_msa_server: bool = True,
    use_potentials: bool = False,
    compute_strain: bool = False,
    do_pose_refinement: bool = False,
    name: str = "Cofolding Workflow",
    model: str = stjames.CofoldingModel.BOLTZ_2.value,
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a protein cofolding workflow to the API.

    :param initial_protein_sequences: The sequences of the proteins to be cofolded.
    :param initial_smiles_list: A list of SMILES strings for the ligands to be cofolded with.
    :param ligand_binding_affinity_index: The index of the ligand for which to compute
        the binding affinity.
    :param use_msa_server: Whether to use the MSA server for the computation.
    :param use_potentials: Whether to use potentials for the computation.
    :param compute_strain: Whether to compute the strain of the pose
        (if `pose_refinement` is enabled).
    :param do_pose_refinement: Whether to optimize non-rotatable bonds in output poses.
    :param name: The name of the workflow.
    :param model: The model to use for the computation.
    :param folder_uuid: The UUID of the folder to store the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    workflow = stjames.ProteinCofoldingWorkflow(
        use_msa_server=use_msa_server,
        use_potentials=use_potentials,
        model=model,
        ligand_binding_affinity_index=ligand_binding_affinity_index,
        initial_smiles_list=initial_smiles_list,
        initial_protein_sequences=initial_protein_sequences,
        do_pose_refinement=do_pose_refinement,
        compute_strain=compute_strain,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "protein_cofolding",
        "workflow_data": workflow.model_dump(mode="json"),
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


__all__ = ["CofoldingResult", "ProteinCofoldingResult", "submit_protein_cofolding_workflow"]
