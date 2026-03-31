"""RBFE graph workflow - build perturbation graphs for relative binding free energy calculations."""

from dataclasses import dataclass
from typing import Literal

import stjames

from ..folder import Folder
from ..molecule import Molecule
from ..types import MoleculeInput
from ..utils import api_client
from .base import Workflow, WorkflowResult, molecule_to_dict, register_result


@dataclass(frozen=True, slots=True)
class RelativeBindingFreeEnergyGraphEdge:
    """An edge in an RBFE perturbation graph."""

    ligand_a: str
    ligand_b: str
    core: list[tuple[int, int]] | None = None
    score: float | None = None
    ddg: float | None = None
    ddg_err: float | None = None
    complex_dg: float | None = None
    complex_dg_err: float | None = None
    solvent_dg: float | None = None
    solvent_dg_err: float | None = None
    vacuum_dg: float | None = None
    vacuum_dg_err: float | None = None
    failed: bool = False
    complex_lambda_values: list[float] | None = None
    complex_overlap_matrix: list[list[float]] | None = None


@register_result("rbfe_graph")
class RelativeBindingFreeEnergyGraphResult(WorkflowResult):
    """Result from an RBFE graph construction workflow."""

    _stjames_class = stjames.RBFEGraphWorkflow

    def __repr__(self) -> str:
        g = self._workflow.graph
        n_edges = len(g.edges) if g else 0
        n_ligands = len(self._workflow.ligands)
        return f"<RelativeBindingFreeEnergyGraphResult edges={n_edges} ligands={n_ligands}>"

    @property
    def ligands(self) -> dict[str, Molecule]:
        """Ligand molecules keyed by identifier."""
        return {k: Molecule.from_stjames(v) for k, v in self._workflow.ligands.items()}

    @property
    def graph(self) -> dict | None:
        """
        The constructed RBFE perturbation graph as a dict, or None if not yet computed.

        Pass directly to ``submit_rbfe_perturbation_workflow(graph=...)``.
        """
        return g.model_dump(mode="json") if (g := self._workflow.graph) else None

    @property
    def edges(self) -> list[RelativeBindingFreeEnergyGraphEdge]:
        """Graph edges. Empty list if graph is not yet built."""
        if not (g := self._workflow.graph):
            return []
        return [
            RelativeBindingFreeEnergyGraphEdge(
                ligand_a=e.ligand_a,
                ligand_b=e.ligand_b,
                core=e.core,
                score=e.score,
                ddg=e.ddg,
                ddg_err=e.ddg_err,
                complex_dg=e.complex_dg,
                complex_dg_err=e.complex_dg_err,
                solvent_dg=e.solvent_dg,
                solvent_dg_err=e.solvent_dg_err,
                vacuum_dg=e.vacuum_dg,
                vacuum_dg_err=e.vacuum_dg_err,
                failed=e.failed,
                complex_lambda_values=e.complex_lambda_values,
                complex_overlap_matrix=e.complex_overlap_matrix,
            )
            for e in g.edges
        ]


def submit_relative_binding_free_energy_graph_workflow(
    ligands: dict[str, MoleculeInput],
    mode: Literal["greedy", "star_map"] = "greedy",
    hub_compound_id: str | None = None,
    greedy_scoring: Literal["best", "jaccard", "dummy_atoms"] = "best",
    greedy_k_min_cut: int = 3,
    refine_cutoff: float | None = None,
    name: str = "RBFE Graph",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
) -> Workflow:
    """
    Submits an RBFE graph construction workflow to the API.

    Builds a perturbation graph connecting ligands for relative binding free
    energy (RBFE) calculations.

    :param ligands: Dictionary mapping ligand identifiers to molecules.
    :param mode: Graph construction strategy: ``"greedy"`` or ``"star_map"``.
    :param hub_compound_id: Ligand identifier for the hub when ``mode="star_map"``.
    :param greedy_scoring: Edge scoring heuristic for greedy mode:
        ``"best"``, ``"jaccard"``, or ``"dummy_atoms"``.
    :param greedy_k_min_cut: Target edge-connectivity for greedy augmentation. Must be > 0.
    :param refine_cutoff: Optional MCS similarity cutoff for graph refinement.
    :param name: Name of the workflow.
    :param folder_uuid: UUID of the folder to place the workflow in.
    :param folder: Folder object to store the workflow in.
    :param max_credits: Maximum number of credits to use for the workflow.
    :param webhook_url: URL that Rowan will POST to when the workflow completes.
    :returns: Workflow object representing the submitted workflow.
    :raises ValueError: If both folder and folder_uuid are provided.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder:
        folder_uuid = folder.uuid

    ligands_dict = {k: molecule_to_dict(v) for k, v in ligands.items()}

    workflow = stjames.RBFEGraphWorkflow(
        ligands=ligands_dict,
        mode=mode,
        hub_compound_id=hub_compound_id,
        greedy_scoring=greedy_scoring,
        greedy_k_min_cut=greedy_k_min_cut,
        refine_cutoff=refine_cutoff,
    )

    data = {
        "workflow_type": "rbfe_graph",
        "workflow_data": workflow.model_dump(mode="json"),
        "name": name,
        "folder_uuid": folder_uuid,
        "max_credits": max_credits,
        "webhook_url": webhook_url,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())
