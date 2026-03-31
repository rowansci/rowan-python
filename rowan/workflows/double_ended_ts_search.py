"""Double-ended TS search workflow - find transition states between reactant and product."""

from dataclasses import dataclass
from typing import Any

import stjames
from stjames.optimization.freezing_string_method import FSMSettings

from ..calculation import Calculation, retrieve_calculation
from ..folder import Folder
from ..molecule import Molecule
from ..types import MoleculeInput
from ..utils import api_client
from .base import Workflow, WorkflowResult, molecule_to_dict, register_result
from .constants import to_relative_kcal


@dataclass(frozen=True, slots=True)
class ReactionPathPoint:
    """
    Point along the reaction path.

    :param distance: Distance along the reaction path.
    :param calculation_uuid: UUID of the calculation at this point.
    :param calculation: Populated by get_path_calculations(), None otherwise.
    """

    distance: float
    calculation_uuid: str | None
    calculation: Calculation | None = None


@register_result("double_ended_ts_search")
class DoubleEndedTSSearchResult(WorkflowResult):
    """Result from a double-ended transition state search workflow."""

    _stjames_class = stjames.DoubleEndedTSSearchWorkflow

    def __post_init__(self) -> None:
        super().__post_init__()
        if self.complete:
            if ts_uuid := getattr(self._workflow, "ts_guess_calculation_uuid", None):
                self._cache["ts_calculation"] = retrieve_calculation(ts_uuid)

    def __repr__(self) -> str:
        ts_uuid = self.ts_guess_calculation_uuid
        n_fwd = len(self.forward_path)
        n_bwd = len(self.backward_path)
        return f"<DoubleEndedTSSearchResult ts_uuid={ts_uuid!r} fwd={n_fwd} bwd={n_bwd}>"

    @property
    def ts_guess_calculation_uuid(self) -> str | None:
        """UUID of the transition state guess calculation."""
        return self._workflow.ts_guess_calculation_uuid

    @property
    def ts_guess_calculation(self) -> Calculation | None:
        """The transition state guess Calculation with full molecule data (lazily fetched)."""
        if "ts_calculation" not in self._cache:
            if ts_uuid := self.ts_guess_calculation_uuid:
                self._cache["ts_calculation"] = retrieve_calculation(ts_uuid)
        return self._cache.get("ts_calculation")

    @property
    def ts_molecule(self) -> Molecule | None:
        """The transition state molecule."""
        calc = self.ts_guess_calculation
        return calc.molecule if calc else None

    @property
    def ts_energy(self) -> float | None:
        """Energy of the transition state (Hartree)."""
        mol = self.ts_molecule
        return mol.energy if mol else None

    def get_path_calculations(self) -> list[ReactionPathPoint]:
        """
        Fetch all path point calculations, sorted by distance.

        Combines forward (reactant → TS) and backward (product → TS) points,
        sorted by distance, each with its Calculation populated.
        Results are cached after the first call.

        .. note::
            Makes one API call per path point on first access.

        :returns: List of ReactionPathPoints with calculations populated, sorted by distance.
        """
        if "path_calculations" not in self._cache:
            all_points = sorted(self.forward_path + self.backward_path, key=lambda p: p.distance)
            self._cache["path_calculations"] = [
                ReactionPathPoint(
                    distance=p.distance,
                    calculation_uuid=p.calculation_uuid,
                    calculation=retrieve_calculation(p.calculation_uuid)
                    if p.calculation_uuid
                    else None,
                )
                for p in all_points
            ]
        return self._cache["path_calculations"]

    def get_path_energies(self, relative: bool = False) -> list[float]:
        """
        Get energies along the reaction path, sorted by distance.

        :param relative: If True, return relative energies in kcal/mol (relative to
            the lowest energy point). If False (default), return absolute energies
            in Hartree.
        :returns: List of energies along the reaction path.
        """
        energies: list[float] = [
            p.calculation.energy
            for p in self.get_path_calculations()
            if p.calculation and p.calculation.energy is not None
        ]
        return to_relative_kcal(energies) if relative else energies

    @property
    def forward_path(self) -> list[ReactionPathPoint]:
        """Reaction path points from reactant to TS."""
        return [
            ReactionPathPoint(distance=d, calculation_uuid=uuid)
            for d, uuid in zip(
                self._workflow.forward_string_distances,
                self._workflow.forward_calculation_uuids,
                strict=False,
            )
        ]

    @property
    def backward_path(self) -> list[ReactionPathPoint]:
        """Reaction path points from product to TS."""
        return [
            ReactionPathPoint(distance=d, calculation_uuid=uuid)
            for d, uuid in zip(
                self._workflow.backward_string_distances,
                self._workflow.backward_calculation_uuids,
                strict=False,
            )
        ]


def submit_double_ended_ts_search_workflow(
    reactant: MoleculeInput,
    product: MoleculeInput,
    calculation_settings: stjames.Settings | dict[str, Any] | None = None,
    search_settings: FSMSettings | dict[str, Any] | None = None,
    optimize_inputs: bool = False,
    optimize_ts: bool = True,
    name: str = "Double-Ended TS Search Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
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
    :param folder_uuid: UUID of the folder to place the workflow in.
    :param folder: Folder object to store the workflow in.
    :param max_credits: Maximum number of credits to use for the workflow.
    :param webhook_url: URL that Rowan will POST to when the workflow completes.
    :returns: Workflow object representing the submitted workflow.
    """
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder:
        folder_uuid = folder.uuid

    reactant_dict = molecule_to_dict(reactant)
    product_dict = molecule_to_dict(product)

    workflow = stjames.DoubleEndedTSSearchWorkflow(
        reactant=reactant_dict,
        product=product_dict,
        calculation_settings=calculation_settings,
        search_settings=search_settings,
        optimize_inputs=optimize_inputs,
        optimize_ts=optimize_ts,
    )
    data = {
        "workflow_type": "double_ended_ts_search",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": reactant_dict,
        "name": name,
        "folder_uuid": folder_uuid,
        "max_credits": max_credits,
        "webhook_url": webhook_url,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())
