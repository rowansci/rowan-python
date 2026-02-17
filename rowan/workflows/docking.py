"""Docking workflow - molecular docking to protein targets."""

from dataclasses import dataclass
from typing import Any

import stjames
from rdkit import Chem

from ..protein import Protein
from ..utils import api_client
from .base import RdkitMol, StJamesMolecule, Workflow, WorkflowResult, register_result


@dataclass(frozen=True, slots=True)
class DockingScore:
    """A docking pose with its score."""

    score: float
    pose: str | None = None
    complex_pdb: str | None = None
    posebusters_valid: bool = False
    strain: float | None = None
    rmsd: float | None = None


@register_result("docking")
class DockingResult(WorkflowResult):
    """Result from a docking workflow."""

    _stjames_class = stjames.DockingWorkflow

    def __repr__(self) -> str:
        scores = self.scores
        best = min((s.score for s in scores), default=None)
        return f"<DockingResult poses={len(scores)} best_score={best}>"

    @property
    def scores(self) -> list[DockingScore]:
        """List of docking scores with poses."""
        return [
            DockingScore(
                score=s.score,
                pose=s.pose,
                complex_pdb=s.complex_pdb,
                posebusters_valid=s.posebusters_valid,
                strain=s.strain,
                rmsd=s.rmsd,
            )
            for s in self._workflow.scores
        ]

    @property
    def conformers(self) -> list[str]:
        """UUIDs of optimized conformers."""
        return list(self._workflow.conformers)


def submit_docking_workflow(
    protein: str | Protein,
    pocket: list[list[float]],
    initial_molecule: dict[str, Any] | StJamesMolecule | RdkitMol,
    executable: str = "vina",
    scoring_function: str = "vinardo",
    exhaustiveness: float = 8,
    do_csearch: bool = False,
    do_optimization: bool = False,
    do_pose_refinement: bool = False,
    name: str = "Docking Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a Docking workflow to the API.

    :param protein: The protein to dock. Can be input as a uuid or a Protein object.
    :param pocket: Binding pocket coordinates [[x,y,z], [x,y,z]].
    :param initial_molecule: The initial molecule to be docked.
    :param executable: Which docking implementation to use.
    :param scoring_function: Which docking scoring function to use.
    :param exhaustiveness: Which exhaustiveness to employ.
    :param do_csearch: Whether to perform a conformational search on the ligand.
    :param do_optimization: Whether to perform an optimization on the ligand.
    :param do_pose_refinement: Whether or not to optimize output poses.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted docking workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if isinstance(initial_molecule, StJamesMolecule):
        initial_molecule = initial_molecule.model_dump(mode="json")
    elif isinstance(initial_molecule, Chem.rdchem.Mol | Chem.rdchem.RWMol):
        initial_molecule = StJamesMolecule.from_rdkit(initial_molecule, cid=0)

    if isinstance(protein, Protein):
        protein = protein.uuid

    docking_settings = {
        "executable": executable,
        "exhaustiveness": exhaustiveness,
        "scoring_function": scoring_function,
    }

    workflow = stjames.DockingWorkflow(
        initial_molecule=initial_molecule,
        protein=protein,
        pocket=pocket,
        do_csearch=do_csearch,
        do_optimization=do_optimization,
        do_pose_refinement=do_pose_refinement,
        docking_settings=docking_settings,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "docking",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": initial_molecule,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


__all__ = ["DockingResult", "DockingScore", "submit_docking_workflow"]
