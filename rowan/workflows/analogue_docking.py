"""Analogue docking workflow - dock analogues using a template ligand."""

from typing import Literal

import stjames

from ..calculation import Calculation, retrieve_calculation
from ..folder import Folder
from ..molecule import Molecule
from ..protein import Protein, retrieve_protein
from ..utils import api_client
from .base import (
    StructureInput,
    Workflow,
    WorkflowResult,
    molecule_to_dict,
    register_result,
    require_coordinates,
)
from .docking import DockingScore


@register_result("analogue_docking")
class AnalogueDockingResult(WorkflowResult):
    """Result from an analogue-docking workflow."""

    _stjames_class = stjames.AnalogueDockingWorkflow

    def __repr__(self) -> str:
        scores = self._workflow.analogue_scores
        n = len(scores)
        # Find best (lowest) score across all analogues
        best_score = float("inf")
        best_smiles = None
        for smiles, score_list in scores.items():
            if score_list:
                min_score = min(s.score for s in score_list)
                if min_score < best_score:
                    best_score = min_score
                    best_smiles = smiles
        return f"<AnalogueDockingResult analogues={n} best=({best_score:.2f}, {best_smiles!r})>"

    def __post_init__(self) -> None:
        """Default `complex_pdb` to None on each analogue's scores, then parse."""
        for scores in (self.workflow_data.get("analogue_scores") or {}).values():
            for score in scores:
                score.setdefault("complex_pdb", None)
        super().__post_init__()

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
                    mmgbsa_score=s.mmgbsa_score,
                )
                for s in scores
            ]
            for smiles, scores in self._workflow.analogue_scores.items()
        }

    def get_pose(self, smiles: str, index: int = 0) -> Calculation:
        """
        Fetch a docked ligand pose as a calculation with 3D coordinates.

        :param smiles: SMILES string of the analogue.
        :param index: Index of the pose (0-based, ordered by score). Default 0 (best).
        :returns: Calculation containing the docked ligand molecule with 3D coordinates.
        :raises KeyError: If the SMILES is not found.
        :raises IndexError: If index is out of range.
        :raises ValueError: If the pose has no UUID.
        """
        scores = self.analogue_scores.get(smiles)
        if scores is None:
            raise KeyError(f"Analogue '{smiles}' not found")
        if index < 0 or index >= len(scores):
            raise IndexError(f"Pose index {index} out of range (0-{len(scores) - 1})")

        uuid = scores[index].pose
        if not uuid:
            raise ValueError(f"Pose {index} for '{smiles}' has no UUID")

        cache_key = f"pose_{smiles}_{index}"
        if cache_key not in self._cache:
            self._cache[cache_key] = retrieve_calculation(uuid)
        return self._cache[cache_key]

    @property
    def best_poses(self) -> dict[str, Molecule]:
        """
        Best docked pose per analogue, keyed by SMILES.

        Fetches the final geometry from each analogue's top-scoring pose.
        Analogues with no successful poses are excluded. When ``analogue_names``
        were supplied at submission, each pose's ``name`` is set to the analogue
        name, so poses can be re-keyed by name (e.g. ``{p.name: p for p in ...}``).

        :returns: Dictionary mapping SMILES to docked ``Molecule`` with 3D coordinates.
        """
        # analogue_names is empty or parallel to analogues; zip drops the tail when empty.
        name_by_smiles = dict(
            zip(self._workflow.analogues, self._workflow.analogue_names, strict=False)
        )
        poses = {}
        for smiles, scores in self.analogue_scores.items():
            if not scores:
                continue
            pose = self.get_pose(smiles).molecules[-1]
            # Keep any name the pose already carries; otherwise fall back to analogue_names.
            if not pose.name:
                pose.name = name_by_smiles.get(smiles)
            poses[smiles] = pose
        return poses

    def get_poses(self, smiles: str) -> list[Calculation]:
        """
        Fetch all docked ligand poses as calculations with 3D coordinates.

        :param smiles: SMILES string of the analogue.
        :returns: List of Calculations for each pose (ordered by score).
        :raises KeyError: If the SMILES is not found.
        """
        scores = self.analogue_scores.get(smiles)
        if scores is None:
            raise KeyError(f"Analogue '{smiles}' not found")

        poses: list[Calculation] = []
        for i, score in enumerate(scores):
            if score.pose:
                poses.append(self.get_pose(smiles, i))
        return poses

    def get_complex(self, smiles: str, index: int = 0) -> Protein:
        """
        Fetch a protein-ligand complex structure for a specific analogue.

        :param smiles: SMILES string of the analogue.
        :param index: Index of the pose (0-based, ordered by score). Default 0 (best).
        :returns: Protein object with the full protein-ligand complex.
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
        :returns: List of Protein objects for each complex (ordered by score).
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
    initial_molecule: StructureInput,
    protein: str | Protein,
    analogue_names: list[str] | None = None,
    scoring_function: Literal["vina", "vinardo"] = "vinardo",
    exhaustiveness: int = 8,
    max_poses: int = 4,
    num_conformers_per_analogue: int = 100,
    require_posebusters: bool = False,
    run_local_optimization: bool = False,
    name: str = "Analogue Docking Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
    is_draft: bool = False,
) -> Workflow:
    """
    Submits an analogue-docking workflow to the API.

    :param analogues: SMILES strings to dock.
    :param initial_molecule: Template to which to align molecules to.
    :param protein: Protein to dock. Can be input as a uuid or a Protein object.
    :param analogue_names: Names parallel to `analogues`. When provided, each docked
        pose's `Molecule.name` is set to its analogue name (see `best_poses`), and must
        be the same length as `analogues`. If omitted, poses are identified by SMILES.
    :param scoring_function: Docking scoring function: "vina" or "vinardo".
    :param exhaustiveness: How many times Vina attempts to find a pose for each conformer.
    :param max_poses: Maximum number of poses generated per input conformer.
    :param num_conformers_per_analogue: Maximum number of conformers to generate per analogue.
    :param require_posebusters: Filter conformers based on PoseBusters validity before docking.
    :param run_local_optimization: Whether to run a local opt in docking pocket or just score.
    :param name: Name of the workflow.
    :param folder_uuid: UUID of the folder to place the workflow in.
    :param folder: Folder object to store the workflow in.
    :param max_credits: Maximum number of credits to use for the workflow.
    :param webhook_url: URL that Rowan will POST to when the workflow completes.
    :param is_draft: If True, submit the workflow as a draft without starting execution.
    :returns: Workflow object representing the submitted analogue-docking workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    require_coordinates(initial_molecule)
    if analogue_names is not None and len(analogue_names) != len(analogues):
        raise ValueError(
            f"analogue_names must match analogues in length "
            f"(got {len(analogue_names)} names for {len(analogues)} analogues)."
        )
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder:
        folder_uuid = folder.uuid
    docking_settings = {
        "executable": "vina",
        "scoring_function": scoring_function,
        "exhaustiveness": exhaustiveness,
        "max_poses": max_poses,
    }

    mol_dict = molecule_to_dict(initial_molecule)

    if isinstance(protein, Protein):
        protein = protein.uuid

    workflow = stjames.AnalogueDockingWorkflow(
        analogues=analogues,
        analogue_names=analogue_names or [],
        initial_molecule=mol_dict,
        protein=protein,
        docking_settings=docking_settings,
        num_conformers_per_analogue=num_conformers_per_analogue,
        require_posebusters=require_posebusters,
        run_local_optimization=run_local_optimization,
    )

    data = {
        "initial_molecule": mol_dict,
        "workflow_type": "analogue_docking",
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
