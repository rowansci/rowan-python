"""Protein binder design workflow - generate protein binders."""

from dataclasses import dataclass
from typing import Any

import stjames

from ..utils import api_client
from .base import Message, Workflow, WorkflowResult, parse_messages, register_result

# Re-export protocol enum from stjames
BinderProtocol = stjames.BoltzGenProtocol


@dataclass(frozen=True, slots=True)
class BinderScores:
    """
    Scores for a generated protein binder design.

    :param iptm: Interface predicted TM-score (0-1, higher is better).
    :param design_ptm: Predicted TM-score for the designed binder (0-1).
    :param quality_score: Overall quality score (0-1, higher is better).
    :param bb_rmsd: Backbone RMSD compared to initial structure (Angstrom).
    :param loop: Fraction of residues in loop conformation.
    :param helix: Fraction of residues in helix conformation.
    :param sheet: Fraction of residues in sheet conformation.
    :param liability_score: Liability score (lower is better).
    :param liability_num_violations: Number of liability violations.
    :param liability_high_severity_violations: Number of high-severity liability violations.
    :param min_interaction_pae: Minimum predicted aligned error at the interface.
    :param delta_sasa_refolded: Change in solvent-accessible surface area upon binding (A^2).
    :param plip_hbonds_refolded: Number of hydrogen bonds at the interface.
    :param plip_saltbridge_refolded: Number of salt bridges at the interface.
    :param num_tokens: Number of tokens in the design.
    :param design_hydrophobicity: Hydrophobicity of the designed binder.
    :param num_filters_passed: Number of quality filters passed.
    """

    iptm: float | None = None
    design_ptm: float | None = None
    quality_score: float | None = None
    bb_rmsd: float | None = None
    loop: float | None = None
    helix: float | None = None
    sheet: float | None = None
    liability_score: float | None = None
    liability_num_violations: int | None = None
    liability_high_severity_violations: int | None = None
    min_interaction_pae: float | None = None
    delta_sasa_refolded: float | None = None
    plip_hbonds_refolded: int | None = None
    plip_saltbridge_refolded: int | None = None
    num_tokens: int | None = None
    design_hydrophobicity: float | None = None
    num_filters_passed: int | None = None


@dataclass(frozen=True, slots=True)
class ProteinBinder:
    """
    Generated protein binder design.

    :param bound_structure_uuid: UUID of the bound structure (binder + target complex).
    :param sequence: Amino acid sequence of the designed binder.
    :param scores: Detailed scores for the binder design.
    """

    bound_structure_uuid: str | None = None
    sequence: str | None = None
    scores: BinderScores | None = None

    # Convenience accessors for common scores
    @property
    def iptm(self) -> float | None:
        """Interface predicted TM-score (0-1, higher is better)."""
        return self.scores.iptm if self.scores else None

    @property
    def quality_score(self) -> float | None:
        """Overall quality score (0-1, higher is better)."""
        return self.scores.quality_score if self.scores else None


@register_result("protein_binder_design")
class ProteinBinderDesignResult(WorkflowResult):
    """Result from a protein-binder-design workflow."""

    _stjames_class = stjames.ProteinBinderDesignWorkflow

    def __repr__(self) -> str:
        binders = self.generated_binders
        n = len(binders)
        if binders:
            best = max(binders, key=lambda b: b.iptm or 0)
            return f"<ProteinBinderDesignResult binders={n} best_iptm={best.iptm}>"
        return f"<ProteinBinderDesignResult binders={n}>"

    @property
    def generated_binders(self) -> list[ProteinBinder]:
        """Generated protein binder designs, sorted by quality score."""
        binders: list[ProteinBinder] = []
        for b in self._workflow.generated_binders:
            scores_dict = getattr(b, "scores", None) or {}
            scores = BinderScores(
                iptm=scores_dict.get("iptm"),
                design_ptm=scores_dict.get("design_ptm"),
                quality_score=scores_dict.get("quality_score"),
                bb_rmsd=scores_dict.get("bb_rmsd"),
                loop=scores_dict.get("loop"),
                helix=scores_dict.get("helix"),
                sheet=scores_dict.get("sheet"),
                liability_score=scores_dict.get("liability_score"),
                liability_num_violations=scores_dict.get("liability_num_violations"),
                liability_high_severity_violations=scores_dict.get(
                    "liability_high_severity_violations"
                ),
                min_interaction_pae=scores_dict.get("min_interaction_pae"),
                delta_sasa_refolded=scores_dict.get("delta_sasa_refolded"),
                plip_hbonds_refolded=scores_dict.get("plip_hbonds_refolded"),
                plip_saltbridge_refolded=scores_dict.get("plip_saltbridge_refolded"),
                num_tokens=scores_dict.get("num_tokens"),
                design_hydrophobicity=scores_dict.get("design_hydrophobicity"),
                num_filters_passed=scores_dict.get("num_filters_passed"),
            )
            binders.append(
                ProteinBinder(
                    bound_structure_uuid=getattr(b, "bound_structure", None),
                    sequence=getattr(b, "sequence", None),
                    scores=scores,
                )
            )
        return binders

    @property
    def messages(self) -> list[Message]:
        """Any messages or warnings from the workflow."""
        return parse_messages(getattr(self._workflow, "messages", None))


def submit_protein_binder_design_workflow(
    binder_design_input: dict[str, Any],
    protocol: BinderProtocol | str = BinderProtocol.PROTEIN_ANYTHING,
    num_designs: int = 10,
    budget: int = 2,
    name: str = "Protein Binder Design Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a protein-binder-design workflow to the API.

    :param binder_design_input: Input specification for the binder design (BoltzGenInput format).
    :param protocol: Design protocol to use. Options:
        - PROTEIN_ANYTHING: Design a protein binder
        - PEPTIDE_ANYTHING: Design a peptide binder
        - PROTEIN_SMALL_MOLECULE: Design a protein that binds a small molecule
        - NANOBODY_ANYTHING: Design a nanobody binder
    :param num_designs: Number of designs to generate.
    :param budget: Number of designs to return in the final diversity-optimized set.
    :param name: Name of the workflow.
    :param folder_uuid: UUID of the folder to place the workflow in.
    :param max_credits: Maximum number of credits to use for the workflow.
    :returns: Workflow object representing the submitted workflow.
    :raises ValueError: If protocol is not a valid BinderProtocol.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if isinstance(protocol, BinderProtocol):
        protocol = protocol.value
    elif isinstance(protocol, str):
        valid = [p.value for p in BinderProtocol]
        if protocol not in valid:
            raise ValueError(f"Invalid protocol '{protocol}'. Must be one of: {', '.join(valid)}")

    workflow = stjames.ProteinBinderDesignWorkflow(
        binder_design_input=binder_design_input,
        binder_design_settings={
            "protocol": protocol,
            "num_designs": num_designs,
            "budget": budget,
        },
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "protein_binder_design",
        "workflow_data": workflow.model_dump(mode="json"),
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())
