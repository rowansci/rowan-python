"""Protein binder design workflow - generate protein binders."""

from dataclasses import dataclass
from typing import Any

import stjames

from ..utils import api_client
from .base import Workflow, WorkflowResult, register_result


@dataclass(frozen=True, slots=True)
class ProteinBinder:
    """A generated protein binder design."""

    sequence: str
    average_plddt: float | None = None
    ptm: float | None = None
    iptm: float | None = None


@register_result("protein_binder_design")
class ProteinBinderDesignResult(WorkflowResult):
    """Result from a protein binder design workflow."""

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
        """Generated protein binder designs."""
        return [
            ProteinBinder(
                sequence=b.sequence,
                average_plddt=getattr(b, "plddt", None),
                ptm=getattr(b, "ptm", None),
                iptm=getattr(b, "iptm", None),
            )
            for b in self._workflow.generated_binders
        ]


def submit_protein_binder_design_workflow(
    binder_design_input: dict[str, Any],
    protocol: str = "protein-anything",
    num_designs: int = 10,
    budget: int = 2,
    name: str = "Protein Binder Design Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a protein binder design workflow to the API.

    :param binder_design_input: Input specification for the binder design (BoltzGenInput format).
    :param protocol: The protocol to use.
    :param num_designs: Number of designs to generate.
    :param budget: Number of designs to return in the final diversity-optimized set.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
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


__all__ = ["ProteinBinder", "ProteinBinderDesignResult", "submit_protein_binder_design_workflow"]
