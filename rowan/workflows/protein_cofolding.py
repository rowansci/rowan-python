"""Protein cofolding workflow - predict protein-protein and protein-ligand complexes."""

from dataclasses import dataclass

import stjames

from ..protein import Protein, retrieve_protein
from ..utils import api_client
from .base import Message, Workflow, WorkflowResult, parse_messages, register_result

# Re-export cofolding model enum from stjames
CofoldingModel = stjames.CofoldingModel


@dataclass(frozen=True, slots=True)
class CofoldingScores:
    """Confidence scores for a cofolding prediction."""

    ptm: float | None = None
    """Predicted TM-score (0-1, higher is better)."""

    iptm: float | None = None
    """Interface predicted TM-score (0-1, higher is better)."""

    avg_lddt: float | None = None
    """Average per-residue LDDT confidence (0-1)."""

    confidence_score: float | None = None
    """Overall confidence score (0-1)."""


@dataclass(frozen=True, slots=True)
class AffinityScore:
    """Predicted binding affinity scores."""

    pred_value: float | None = None
    """Predicted binding affinity (ensemble average)."""

    pred_value1: float | None = None
    """Predicted binding affinity (model 1)."""

    pred_value2: float | None = None
    """Predicted binding affinity (model 2)."""

    probability_binary: float | None = None
    """Probability of binding (ensemble average, 0-1)."""

    probability_binary1: float | None = None
    """Probability of binding (model 1, 0-1)."""

    probability_binary2: float | None = None
    """Probability of binding (model 2, 0-1)."""


@dataclass(frozen=True, slots=True)
class CofoldingResult:
    """A single cofolding prediction result."""

    scores: CofoldingScores | None = None
    """Confidence scores for the prediction."""

    affinity_score: AffinityScore | None = None
    """Predicted binding affinity (if computed)."""

    strain: float | None = None
    """Ligand strain energy (if computed)."""

    posebusters_valid: bool | None = None
    """Whether the pose passes PoseBusters validation."""

    lddt: list[float] | None = None
    """Per-residue LDDT confidence scores."""

    pose_uuid: str | None = None
    """UUID of the pose."""

    predicted_structure_uuid: str | None = None
    """UUID of the predicted structure."""

    predicted_refined_structure_uuid: str | None = None
    """UUID of the refined structure (if pose refinement was enabled)."""


@register_result("protein_cofolding")
class ProteinCofoldingResult(WorkflowResult):
    """Result from a protein cofolding workflow."""

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
        return CofoldingScores(
            ptm=s.get("ptm") if isinstance(s, dict) else getattr(s, "ptm", None),
            iptm=s.get("iptm") if isinstance(s, dict) else getattr(s, "iptm", None),
            avg_lddt=s.get("avg_lddt") if isinstance(s, dict) else getattr(s, "avg_lddt", None),
            confidence_score=(
                s.get("confidence_score")
                if isinstance(s, dict)
                else getattr(s, "confidence_score", None)
            ),
        )

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
        """Fetch the predicted structure as a Protein object."""
        uuid = self.predicted_structure_uuid
        if not uuid:
            return None
        if "predicted_structure" not in self._cache:
            self._cache["predicted_structure"] = retrieve_protein(uuid)
        return self._cache["predicted_structure"]

    def get_refined_structure(self) -> Protein | None:
        """Fetch the refined structure as a Protein object (if available)."""
        uuid = self.predicted_refined_structure_uuid
        if not uuid:
            return None
        if "refined_structure" not in self._cache:
            self._cache["refined_structure"] = retrieve_protein(uuid)
        return self._cache["refined_structure"]

    @property
    def predictions(self) -> list[CofoldingResult]:
        """All cofolding predictions."""
        results: list[CofoldingResult] = []
        for r in self._workflow.cofolding_results or []:
            scores_data = getattr(r, "scores", None)
            scores = None
            if scores_data:
                scores = CofoldingScores(
                    ptm=(
                        scores_data.get("ptm")
                        if isinstance(scores_data, dict)
                        else getattr(scores_data, "ptm", None)
                    ),
                    iptm=(
                        scores_data.get("iptm")
                        if isinstance(scores_data, dict)
                        else getattr(scores_data, "iptm", None)
                    ),
                    avg_lddt=(
                        scores_data.get("avg_lddt")
                        if isinstance(scores_data, dict)
                        else getattr(scores_data, "avg_lddt", None)
                    ),
                    confidence_score=(
                        scores_data.get("confidence_score")
                        if isinstance(scores_data, dict)
                        else getattr(scores_data, "confidence_score", None)
                    ),
                )

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
        return parse_messages(getattr(self._workflow, "messages", None))

    def _make_affinity_score(self, a: dict | object) -> AffinityScore:
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
            pred_value=getattr(a, "pred_value", None),
            pred_value1=getattr(a, "pred_value1", None),
            pred_value2=getattr(a, "pred_value2", None),
            probability_binary=getattr(a, "probability_binary", None),
            probability_binary1=getattr(a, "probability_binary1", None),
            probability_binary2=getattr(a, "probability_binary2", None),
        )


def submit_protein_cofolding_workflow(
    initial_protein_sequences: list[str] | None = None,
    initial_dna_sequences: list[str] | None = None,
    initial_rna_sequences: list[str] | None = None,
    initial_smiles_list: list[str] | None = None,
    ligand_binding_affinity_index: int | None = None,
    use_msa_server: bool = True,
    use_potentials: bool = False,
    compute_strain: bool = False,
    do_pose_refinement: bool = False,
    name: str = "Protein–Ligand Co-Folding",
    model: CofoldingModel | str = CofoldingModel.BOLTZ_2,
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a protein cofolding workflow to the API.

    Predicts the 3D structure of protein-protein, protein-ligand, protein-DNA,
    protein-RNA, or other biomolecular complexes.

    :param initial_protein_sequences: Protein sequences to be cofolded.
    :param initial_dna_sequences: DNA sequences to be cofolded.
    :param initial_rna_sequences: RNA sequences to be cofolded.
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
    :raises ValueError: If no protein, DNA, or RNA sequences are provided.
    :raises requests.HTTPError: if the request to the API fails.
    """
    # At least one sequence type is required
    has_protein = initial_protein_sequences and len(initial_protein_sequences) > 0
    has_dna = initial_dna_sequences and len(initial_dna_sequences) > 0
    has_rna = initial_rna_sequences and len(initial_rna_sequences) > 0

    if not (has_protein or has_dna or has_rna):
        raise ValueError(
            "At least one protein, DNA, or RNA sequence is required for cofolding"
        )

    # Validate and convert model
    if isinstance(model, CofoldingModel):
        model_str = model.value
    else:
        valid_models = [m.value for m in CofoldingModel]
        if model not in valid_models:
            raise ValueError(
                f"Invalid model '{model}'. Must be one of: {', '.join(valid_models)}"
            )
        model_str = model

    workflow = stjames.ProteinCofoldingWorkflow(
        use_msa_server=use_msa_server,
        use_potentials=use_potentials,
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


__all__ = [
    "AffinityScore",
    "CofoldingModel",
    "CofoldingResult",
    "CofoldingScores",
    "ProteinCofoldingResult",
    "submit_protein_cofolding_workflow",
]
