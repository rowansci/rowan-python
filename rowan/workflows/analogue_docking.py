"""Analogue docking workflow - dock analogues using a template ligand."""

from typing import Any

import stjames
from rdkit import Chem

from ..protein import Protein, retrieve_protein
from ..utils import api_client
from .base import RdkitMol, StJamesMolecule, Workflow, WorkflowResult, register_result
from .docking import DockingScore


@register_result("analogue_docking")
class AnalogueDockingResult(WorkflowResult):
    """Result from an analogue docking workflow."""

    _stjames_class = stjames.AnalogueDockingWorkflow

    def __repr__(self) -> str:
        scores = self._workflow.analogue_scores
        n = len(scores)
        # Find best (lowest) score across all analogues
        best_score = None
        best_smiles = None
        for smiles, score_list in scores.items():
            if score_list:
                min_score = min(s.score for s in score_list)
                if best_score is None or min_score < best_score:
                    best_score = min_score
                    best_smiles = smiles
        return f"<AnalogueDockingResult analogues={n} best=({best_score}, {best_smiles!r})>"

    @property
    def analogue_scores(self) -> dict[str, list[DockingScore]]:
        """Docking scores for each analogue SMILES."""
        return {
            smiles: [
                DockingScore(
                    score=s.score,
                    pose=s.pose,
                    complex_pdb=s.complex_pdb,
                    posebusters_valid=s.posebusters_valid,
                    strain=s.strain,
                    rmsd=s.rmsd,
                )
                for s in scores
            ]
            for smiles, scores in self._workflow.analogue_scores.items()
        }

    def get_pose(self, smiles: str, index: int = 0) -> Protein:
        """
        Fetch a docked pose structure for a specific analogue.

        :param smiles: SMILES string of the analogue.
        :param index: Index of the pose (0-based, ordered by score). Default 0 (best).
        :return: Protein object with the docked ligand pose.
        :raises KeyError: If the SMILES is not found.
        :raises IndexError: If index is out of range.
        :raises ValueError: If the pose has no structure UUID.
        """
        scores = self.analogue_scores.get(smiles)
        if scores is None:
            raise KeyError(f"Analogue '{smiles}' not found")
        if index < 0 or index >= len(scores):
            raise IndexError(f"Pose index {index} out of range (0-{len(scores) - 1})")

        uuid = scores[index].pose
        if not uuid:
            raise ValueError(f"Pose {index} for '{smiles}' has no structure UUID")

        cache_key = f"pose_{smiles}_{index}"
        if cache_key not in self._cache:
            self._cache[cache_key] = retrieve_protein(uuid)
        return self._cache[cache_key]

    def get_poses(self, smiles: str) -> list[Protein]:
        """
        Fetch all docked pose structures for a specific analogue.

        :param smiles: SMILES string of the analogue.
        :return: List of Protein objects for each pose (ordered by score).
        :raises KeyError: If the SMILES is not found.
        """
        scores = self.analogue_scores.get(smiles)
        if scores is None:
            raise KeyError(f"Analogue '{smiles}' not found")

        poses: list[Protein] = []
        for i, score in enumerate(scores):
            if score.pose:
                poses.append(self.get_pose(smiles, i))
        return poses

    def get_complex(self, smiles: str, index: int = 0) -> Protein:
        """
        Fetch a protein-ligand complex structure for a specific analogue.

        :param smiles: SMILES string of the analogue.
        :param index: Index of the pose (0-based, ordered by score). Default 0 (best).
        :return: Protein object with the full protein-ligand complex.
        :raises KeyError: If the SMILES is not found.
        :raises IndexError: If index is out of range.
        :raises ValueError: If the complex has no structure UUID.
        """
        scores = self.analogue_scores.get(smiles)
        if scores is None:
            raise KeyError(f"Analogue '{smiles}' not found")
        if index < 0 or index >= len(scores):
            raise IndexError(f"Pose index {index} out of range (0-{len(scores) - 1})")

        uuid = scores[index].complex_pdb
        if not uuid:
            raise ValueError(f"Pose {index} for '{smiles}' has no complex UUID")

        cache_key = f"complex_{smiles}_{index}"
        if cache_key not in self._cache:
            self._cache[cache_key] = retrieve_protein(uuid)
        return self._cache[cache_key]

    def get_complexes(self, smiles: str) -> list[Protein]:
        """
        Fetch all protein-ligand complex structures for a specific analogue.

        :param smiles: SMILES string of the analogue.
        :return: List of Protein objects for each complex (ordered by score).
        :raises KeyError: If the SMILES is not found.
        """
        scores = self.analogue_scores.get(smiles)
        if scores is None:
            raise KeyError(f"Analogue '{smiles}' not found")

        complexes: list[Protein] = []
        for i, score in enumerate(scores):
            if score.complex_pdb:
                complexes.append(self.get_complex(smiles, i))
        return complexes


def submit_analogue_docking_workflow(
    analogues: list[str],
    initial_molecule: dict[str, Any] | stjames.Molecule | RdkitMol,
    protein: str | Protein,
    executable: str = "vina",
    scoring_function: str = "vinardo",
    exhaustiveness: float = 8,
    name: str = "Analogue Docking Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits an analogue docking workflow to the API.

    :param analogues: The SMILES strings to dock.
    :param initial_molecule: The template to which to align molecules to.
    :param protein: The protein to dock. Can be input as a uuid or a Protein object.
    :param executable: Which docking implementation to use.
    :param scoring_function: Which docking scoring function to use.
    :param exhaustiveness: Which exhaustiveness to employ.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted batch docking workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    docking_settings = {
        "executable": executable,
        "exhaustiveness": exhaustiveness,
        "scoring_function": scoring_function,
    }

    if isinstance(initial_molecule, StJamesMolecule):
        initial_molecule = initial_molecule.model_dump(mode="json")
    elif isinstance(initial_molecule, Chem.rdchem.Mol | Chem.rdchem.RWMol):
        initial_molecule = stjames.Molecule.from_rdkit(initial_molecule, cid=0).model_dump(
            mode="json"
        )

    if isinstance(protein, Protein):
        protein = protein.uuid

    workflow = stjames.AnalogueDockingWorkflow(
        analogues=analogues,
        initial_molecule=initial_molecule,
        protein=protein,
        docking_settings=docking_settings,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "initial_molecule": initial_molecule,
        "workflow_type": "analogue_docking",
        "workflow_data": workflow.model_dump(serialize_as_any=True, mode="json"),
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


__all__ = ["AnalogueDockingResult", "DockingScore", "submit_analogue_docking_workflow"]
