"""Docking workflow - molecular docking to protein targets."""

from dataclasses import dataclass

import stjames

from ..folder import Folder
from ..protein import Protein, retrieve_protein
from ..types import MoleculeInput
from ..utils import api_client
from .base import Workflow, WorkflowResult, molecule_to_dict, register_result


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
        return self._workflow.conformers

    def get_pose(self, index: int = 0) -> Protein:
        """
        Fetch a docked pose structure.

        :param index: Index of the pose (0-based, ordered by score). Default 0 (best).
        :returns: Protein object with the docked ligand pose.
        :raises IndexError: If index is out of range.
        :raises ValueError: If the pose has no structure UUID.
        """
        scores = self.scores
        if index < 0 or index >= len(scores):
            raise IndexError(f"Pose index {index} out of range (0-{len(scores) - 1})")

        uuid = scores[index].pose
        if not uuid:
            raise ValueError(f"Pose {index} has no structure UUID")

        cache_key = f"pose_{index}"
        if cache_key not in self._cache:
            self._cache[cache_key] = retrieve_protein(uuid)
        return self._cache[cache_key]

    def get_poses(self) -> list[Protein]:
        """
        Fetch all docked pose structures.

        :returns: List of Protein objects for each pose (ordered by score).
        """
        poses: list[Protein] = []
        for i, score in enumerate(self.scores):
            if score.pose:
                poses.append(self.get_pose(i))
        return poses

    def get_complex(self, index: int = 0) -> Protein:
        """
        Fetch a protein-ligand complex structure.

        :param index: Index of the pose (0-based, ordered by score). Default 0 (best).
        :returns: Protein object with the full protein-ligand complex.
        :raises IndexError: If index is out of range.
        :raises ValueError: If the complex has no structure UUID.
        """
        scores = self.scores
        if index < 0 or index >= len(scores):
            raise IndexError(f"Pose index {index} out of range (0-{len(scores) - 1})")

        uuid = scores[index].complex_pdb
        if not uuid:
            raise ValueError(f"Pose {index} has no complex UUID")

        cache_key = f"complex_{index}"
        if cache_key not in self._cache:
            self._cache[cache_key] = retrieve_protein(uuid)
        return self._cache[cache_key]

    def get_complexes(self) -> list[Protein]:
        """
        Fetch all protein-ligand complex structures.

        :returns: List of Protein objects for each complex (ordered by score).
        """
        complexes: list[Protein] = []
        for i, score in enumerate(self.scores):
            if score.complex_pdb:
                complexes.append(self.get_complex(i))
        return complexes


def submit_docking_workflow(
    protein: str | Protein,
    pocket: list[list[float]],
    initial_molecule: MoleculeInput,
    executable: str = "vina",
    scoring_function: str = "vinardo",
    exhaustiveness: float = 8,
    do_csearch: bool = False,
    do_optimization: bool = False,
    do_pose_refinement: bool = True,
    name: str = "Docking Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a docking workflow to the API.

    :param protein: Protein to dock. Can be input as a uuid or a Protein object.
    :param pocket: Binding pocket coordinates [[x,y,z], [x,y,z]].
    :param initial_molecule: Initial molecule to be docked.
    :param executable: Which docking implementation to use.
    :param scoring_function: Which docking scoring function to use.
    :param exhaustiveness: Which exhaustiveness to employ.
    :param do_csearch: Whether to perform a conformational search on the ligand.
    :param do_optimization: Whether to perform an optimization on the ligand.
    :param do_pose_refinement: Whether or not to optimize output poses.
    :param name: Name of the workflow.
    :param folder_uuid: UUID of the folder to place the workflow in.
    :param folder: Folder object to store the workflow in.
    :param max_credits: Maximum number of credits to use for the workflow.
    :returns: Workflow object representing the submitted docking workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder:
        folder_uuid = folder.uuid
    initial_molecule = molecule_to_dict(initial_molecule)

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
