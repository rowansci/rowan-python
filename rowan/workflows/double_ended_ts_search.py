"""Double-ended TS search workflow - find transition states between reactant and product."""

from typing import Any

import stjames
from stjames.optimization.freezing_string_method import FSMSettings

from ..utils import api_client
from .base import StJamesMolecule, Workflow, WorkflowResult, register_result


@register_result("double_ended_ts_search")
class DoubleEndedTSSearchResult(WorkflowResult):
    """Result from a double-ended transition state search workflow."""

    _stjames_class = stjames.DoubleEndedTSSearchWorkflow

    def __repr__(self) -> str:
        ts_uuid = self.ts_guess_calculation_uuid
        n_pts = len(self.distances) if self.distances else 0
        return f"<DoubleEndedTSSearchResult ts_uuid={ts_uuid!r} path_points={n_pts}>"

    @property
    def ts_guess_calculation_uuid(self) -> str | None:
        """UUID of the transition state guess calculation."""
        return self._workflow.ts_guess_calculation_uuid

    @property
    def distances(self) -> list[float] | None:
        """Path distances from reactant to product."""
        return list(self._workflow.distances) if self._workflow.distances else None


def submit_double_ended_ts_search_workflow(
    reactant: dict[str, Any] | StJamesMolecule,
    product: dict[str, Any] | StJamesMolecule,
    calculation_settings: stjames.Settings | dict[str, Any] | None = None,
    search_settings: FSMSettings | dict[str, Any] | None = None,
    optimize_inputs: bool = False,
    optimize_ts: bool = True,
    name: str = "Double-Ended TS Search Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a double-ended transition state search workflow to the API.

    :param reactant: reactant Molecule.
    :param product: product Molecule.
    :param calculation_settings: Settings to use for calculations.
    :param search_settings: settings to use for the transition state search.
    :param optimize_inputs: Whether to optimize the reactant and product before the search.
    :param optimize_ts: Whether to optimize the found transition state.
    :param name: name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: Workflow object representing the submitted workflow.
    """
    workflow = stjames.DoubleEndedTSSearchWorkflow(
        reactant=reactant,
        product=product,
        calculation_settings=calculation_settings,
        search_settings=search_settings,
        optimize_inputs=optimize_inputs,
        optimize_ts=optimize_ts,
    )
    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "double_ended_ts_search",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": reactant
        if isinstance(reactant, dict)
        else reactant.model_dump(mode="json"),
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


__all__ = ["DoubleEndedTSSearchResult", "submit_double_ended_ts_search_workflow"]
