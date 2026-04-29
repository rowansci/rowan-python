"""Protein cofolding workflow - predict protein-protein and protein-ligand complexes."""

from dataclasses import dataclass

import stjames

from ..folder import Folder
from ..protein import Protein, retrieve_protein
from ..utils import api_client
from .base import Message, Workflow, WorkflowResult, parse_messages, register_result

CofoldingModel = stjames.CofoldingModel
# `ConstraintTarget` aliases stjames's `Token` - one position in an input
# (a protein/nucleic-acid residue or a ligand atom) addressable by a constraint.
ConstraintTarget = stjames.workflows.protein_cofolding.Token
ContactConstraint = stjames.workflows.protein_cofolding.ContactConstraint
PocketConstraint = stjames.workflows.protein_cofolding.PocketConstraint


@dataclass(frozen=True, slots=True)
class CofoldingScores:
    """
    Confidence scores for a cofolding prediction.

    :param ptm: Predicted TM-score (0-1, higher is better).
    :param iptm: Interface predicted TM-score (0-1, higher is better).
    :param avg_lddt: Average per-residue LDDT confidence (0-1).
    :param confidence_score: Overall confidence score (0-1).
    """

    ptm: float | None = None
    iptm: float | None = None
    avg_lddt: float | None = None
    confidence_score: float | None = None


@dataclass(frozen=True, slots=True)
class AffinityScore:
    """
    Predicted binding affinity scores.

    :param pred_value: Predicted binding affinity (ensemble average).
    :param pred_value1: Predicted binding affinity (model 1).
    :param pred_value2: Predicted binding affinity (model 2).
    :param probability_binary: Probability of binding (ensemble average, 0-1).
    :param probability_binary1: Probability of binding (model 1, 0-1).
    :param probability_binary2: Probability of binding (model 2, 0-1).
    """

    pred_value: float | None = None
    pred_value1: float | None = None
    pred_value2: float | None = None
    probability_binary: float | None = None
    probability_binary1: float | None = None
    probability_binary2: float | None = None


@dataclass(frozen=True, slots=True)
class CofoldingResult:
    """
    Single cofolding prediction result.

    :param scores: Confidence scores for the prediction.
    :param affinity_score: Predicted binding affinity (if computed).
    :param strain: Ligand strain energy (if computed).
    :param posebusters_valid: Whether the pose passes PoseBusters validation.
    :param lddt: Per-residue LDDT confidence scores.
    :param pose_uuid: UUID of the pose.
    :param predicted_structure_uuid: UUID of the predicted structure.
    :param predicted_refined_structure_uuid: UUID of the refined structure (if refinement was run).
    """

    scores: CofoldingScores | None = None
    affinity_score: AffinityScore | None = None
    strain: float | None = None
    posebusters_valid: bool | None = None
    lddt: list[float] | None = None
    pose_uuid: str | None = None
    predicted_structure_uuid: str | None = None
    predicted_refined_structure_uuid: str | None = None


@register_result("protein_cofolding")
class ProteinCofoldingResult(WorkflowResult):
    """Result from a protein-cofolding workflow."""

    _stjames_class = stjames.ProteinCofoldingWorkflow

    def __repr__(self) -> str:
        n = len(self.predictions)
        scores = self.scores
        iptm = scores.iptm if scores else None
        return f"<ProteinCofoldingResult predictions={n} iptm={iptm}>"

    # --- Top-level results (best/primary prediction) ---

    @property
    def scores(self) -> CofoldingScores | None:
        """Confidence scores for the primary prediction."""
        s = getattr(self._workflow, "scores", None)
        if not s:
            return None
        return self._make_cofolding_scores(s)

    @property
    def affinity_score(self) -> AffinityScore | None:
        """Predicted binding affinity for the primary prediction."""
        a = getattr(self._workflow, "affinity_score", None)
        if not a:
            return None
        return self._make_affinity_score(a)

    @property
    def strain(self) -> float | None:
        """Ligand strain energy for the primary prediction."""
        return getattr(self._workflow, "strain", None)

    @property
    def posebusters_valid(self) -> bool | None:
        """Whether the primary pose passes PoseBusters validation."""
        return getattr(self._workflow, "posebusters_valid", None)

    @property
    def lddt(self) -> list[float] | None:
        """Per-residue LDDT confidence scores for the primary prediction."""
        lddt = getattr(self._workflow, "lddt", None)
        return list(lddt) if lddt else None

    @property
    def predicted_structure_uuid(self) -> str | None:
        """UUID of the predicted structure."""
        return getattr(self._workflow, "predicted_structure_uuid", None)

    @property
    def predicted_refined_structure_uuid(self) -> str | None:
        """UUID of the refined structure (if pose refinement was enabled)."""
        return getattr(self._workflow, "predicted_refined_structure_uuid", None)

    def get_predicted_structure(self) -> Protein | None:
        """Fetch the predicted structure as a Protein object.

        .. note::
            Makes one API call on first access.
            Results are cached. Call clear_cache() to refresh.
        """
        if not (uuid := self.predicted_structure_uuid):
            return None
        if "predicted_structure" not in self._cache:
            self._cache["predicted_structure"] = retrieve_protein(uuid)
        return self._cache["predicted_structure"]

    def get_refined_structure(self) -> Protein | None:
        """Fetch the refined structure as a Protein object (if available).

        .. note::
            Makes one API call on first access.
            Results are cached. Call clear_cache() to refresh.
        """
        if not (uuid := self.predicted_refined_structure_uuid):
            return None
        if "refined_structure" not in self._cache:
            self._cache["refined_structure"] = retrieve_protein(uuid)
        return self._cache["refined_structure"]

    @property
    def predictions(self) -> list[CofoldingResult]:
        """All cofolding predictions."""
        results: list[CofoldingResult] = []
        for r in getattr(self._workflow, "cofolding_results", None) or []:
            scores_data = getattr(r, "scores", None)
            scores = self._make_cofolding_scores(scores_data) if scores_data else None

            affinity = getattr(r, "affinity_score", None)
            affinity_score = self._make_affinity_score(affinity) if affinity else None

            lddt = getattr(r, "lddt", None)

            results.append(
                CofoldingResult(
                    scores=scores,
                    affinity_score=affinity_score,
                    strain=getattr(r, "strain", None),
                    posebusters_valid=getattr(r, "posebusters_valid", None),
                    lddt=list(lddt) if lddt else None,
                    pose_uuid=getattr(r, "pose", None),
                    predicted_structure_uuid=getattr(r, "predicted_structure_uuid", None),
                    predicted_refined_structure_uuid=getattr(
                        r, "predicted_refined_structure_uuid", None
                    ),
                )
            )
        return results

    @property
    def cofolding_results(self) -> list[CofoldingResult]:
        """Alias for predictions (matches API response field name)."""
        return self.predictions

    @property
    def messages(self) -> list[Message]:
        """Any messages or warnings from the workflow (e.g., stereochemistry issues)."""
        return parse_messages(self._workflow.messages)

    @staticmethod
    def _make_cofolding_scores(s: dict | stjames.CofoldingScores) -> CofoldingScores:
        """Convert cofolding scores data to CofoldingScores dataclass."""
        if isinstance(s, dict):
            return CofoldingScores(
                ptm=s.get("ptm"),
                iptm=s.get("iptm"),
                avg_lddt=s.get("avg_lddt"),
                confidence_score=s.get("confidence_score"),
            )
        return CofoldingScores(
            ptm=s.ptm,
            iptm=s.iptm,
            avg_lddt=s.avg_lddt,
            confidence_score=s.confidence_score,
        )

    @staticmethod
    def _make_affinity_score(a: dict | stjames.AffinityScore) -> AffinityScore:
        """Convert affinity score data to AffinityScore dataclass."""
        if isinstance(a, dict):
            return AffinityScore(
                pred_value=a.get("pred_value"),
                pred_value1=a.get("pred_value1"),
                pred_value2=a.get("pred_value2"),
                probability_binary=a.get("probability_binary"),
                probability_binary1=a.get("probability_binary1"),
                probability_binary2=a.get("probability_binary2"),
            )
        return AffinityScore(
            pred_value=a.pred_value,
            pred_value1=a.pred_value1,
            pred_value2=a.pred_value2,
            probability_binary=a.probability_binary,
            probability_binary1=a.probability_binary1,
            probability_binary2=a.probability_binary2,
        )


def submit_protein_cofolding_workflow(
    initial_protein_sequences: list[str] | None = None,
    initial_dna_sequences: list[str] | None = None,
    initial_rna_sequences: list[str] | None = None,
    initial_smiles_list: list[str] | None = None,
    ligand_binding_affinity_index: int | None = None,
    use_msa_server: bool = True,
    use_potentials: bool = False,
    contact_constraints: list[ContactConstraint] | None = None,
    pocket_constraints: list[PocketConstraint] | None = None,
    num_samples: int | None = None,
    compute_strain: bool = False,
    do_pose_refinement: bool = False,
    name: str = "Protein-Ligand Co-Folding",
    model: CofoldingModel | str = CofoldingModel.BOLTZ_2,
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
    is_draft: bool = False,
) -> Workflow:
    """
    Submits a protein-cofolding workflow to the API.

    Predicts the 3D structure of protein-protein, protein-ligand, protein-DNA,
    protein-RNA, or other biomolecular complexes.

    See `examples/protein_cofolding_with_constraints.py` for a worked example
    of using `ConstraintTarget`, `ContactConstraint`, and `PocketConstraint`
    (Boltz-2 only).

    :param initial_protein_sequences: Protein sequences to be cofolded.
    :param initial_dna_sequences: DNA sequences to be cofolded.
    :param initial_rna_sequences: RNA sequences to be cofolded.
    :param initial_smiles_list: List of SMILES strings for the ligands to be cofolded with.
    :param ligand_binding_affinity_index: Index of the ligand for which to compute
        the binding affinity.
    :param use_msa_server: Whether to use the MSA server for the computation.
    :param use_potentials: Whether to use potentials (inference-time steering) with Boltz.
    :param contact_constraints: Boltz contact constraints between two tokens.
    :param pocket_constraints: Boltz pocket constraints between a binder and contact tokens.
    :param num_samples: Number of diffusion samples to generate. If None, uses the model default.
    :param compute_strain: Whether to compute the strain of the pose
        (if `pose_refinement` is enabled).
    :param do_pose_refinement: Whether to optimize non-rotatable bonds in output poses.
    :param name: Name of the workflow.
    :param model: Model to use for the computation.
    :param folder_uuid: UUID of the folder to store the workflow in.
    :param folder: Folder object to store the workflow in.
    :param max_credits: Maximum number of credits to use for the workflow.
    :param webhook_url: URL that Rowan will POST to when the workflow completes.
    :param is_draft: If True, submit the workflow as a draft without starting execution.
    :returns: Workflow object representing the submitted workflow.
    :raises ValueError: If no protein, DNA, or RNA sequences are provided.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder:
        folder_uuid = folder.uuid
    # At least one sequence type is required
    if not (initial_protein_sequences or initial_dna_sequences or initial_rna_sequences):
        raise ValueError("At least one protein, DNA, or RNA sequence is required for cofolding")

    # Validate and convert model
    if isinstance(model, CofoldingModel):
        model_str = model.value
    else:
        valid_models = [m.value for m in CofoldingModel]
        if model not in valid_models:
            raise ValueError(f"Invalid model '{model}'. Must be one of: {', '.join(valid_models)}")
        model_str = model

    workflow = stjames.ProteinCofoldingWorkflow(
        use_msa_server=use_msa_server,
        use_potentials=use_potentials,
        contact_constraints=contact_constraints or [],
        pocket_constraints=pocket_constraints or [],
        num_samples=num_samples,
        model=model_str,
        ligand_binding_affinity_index=ligand_binding_affinity_index,
        initial_smiles_list=initial_smiles_list,
        initial_protein_sequences=initial_protein_sequences or [],
        initial_dna_sequences=initial_dna_sequences or [],
        initial_rna_sequences=initial_rna_sequences or [],
        do_pose_refinement=do_pose_refinement,
        compute_strain=compute_strain,
    )

    data = {
        "workflow_type": "protein_cofolding",
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
