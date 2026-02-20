"""Analogue docking workflow - dock analogues using a template ligand."""

from typing import Any

import stjames
from rdkit import Chem

from ..protein import Protein
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
