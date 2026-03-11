"""Scan workflow - perform potential energy surface scans."""

from typing import Any

import stjames

from ..calculation import Calculation, retrieve_calculation
from ..utils import api_client
from .base import (
    Message,
    MoleculeInput,
    Workflow,
    WorkflowResult,
    molecule_to_dict,
    parse_messages,
    register_result,
)
from .constants import HARTREE_TO_KCAL


@register_result("scan")
class ScanResult(WorkflowResult):
    """Result from a scan workflow."""

    _stjames_class = stjames.ScanWorkflow

    def __repr__(self) -> str:
        n = len(self.scan_point_uuids)
        return f"<ScanResult points={n}>"

    @property
    def scan_point_uuids(self) -> list[str]:
        """UUIDs of scan point calculations."""
        return [uuid for uuid in self._workflow.scan_points if uuid]

    @property
    def scan_points(self) -> list[Calculation]:
        """All scan point calculations.

        Note: Makes one API call per scan point on first access.
        Results are cached. Call clear_cache() to refresh.
        """
        if "scan_points" not in self._cache:
            self._cache["scan_points"] = [
                retrieve_calculation(uuid) for uuid in self.scan_point_uuids
            ]
        return self._cache["scan_points"]

    def get_energies(self, relative: bool = False) -> list[tuple[float, float | None]]:
        """
        Get scan coordinate values paired with energies.

        :param relative: If True, return relative energies in kcal/mol (relative to
            the lowest energy point). If False (default), return absolute energies
            in Hartree.
        :returns: List of (coordinate, energy) tuples. Coordinate is the scanned
            value (e.g., bond distance in Angstrom, angle in degrees).
        """
        scan_settings = self._workflow.scan_settings
        calculations = self.scan_points

        # Compute coordinate values
        if not scan_settings or len(scan_settings) == 0:
            coords = [float(i) for i in range(len(calculations))]
        else:
            settings = scan_settings[0]
            start = settings.start
            stop = settings.stop
            num = settings.num
            coords = [start + i * (stop - start) / (num - 1) for i in range(num)]

        # Get energies
        energies: list[float | None] = [calc.energy for calc in calculations]

        if relative:
            valid_energies = [e for e in energies if e is not None]
            if valid_energies:
                min_energy = min(valid_energies)
                energies = [
                    (e - min_energy) * HARTREE_TO_KCAL if e is not None else None for e in energies
                ]

        return list(zip(coords, energies, strict=True))

    @property
    def messages(self) -> list[Message]:
        """Any messages or warnings from the workflow."""
        return parse_messages(getattr(self._workflow, "messages", None))


def submit_scan_workflow(
    initial_molecule: MoleculeInput,
    scan_settings: stjames.ScanSettings | dict[str, Any] | None = None,
    calculation_engine: str | None = None,
    calculation_method: stjames.Method | str = "uma_m_omol",
    wavefront_propagation: bool = True,
    name: str = "Scan Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a scan workflow to the API.

    :param initial_molecule: Molecule to scan.
    :param scan_settings: Scan settings.
    :param calculation_engine: Engine to use for the calculation.
    :param calculation_method: Method to use for the calculation.
    :param wavefront_propagation: Whether to use wavefront propagation in the scan.
    :param name: Name of the workflow.
    :param folder_uuid: UUID of the folder to store the workflow in.
    :param max_credits: Maximum number of credits to use for the workflow.
    :returns: Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    initial_molecule = molecule_to_dict(initial_molecule)

    if isinstance(calculation_method, str):
        calculation_method = stjames.Method(calculation_method)

    calc_settings = stjames.Settings(
        method=calculation_method,
        tasks=["optimize"],
        corrections=[],
        mode="auto",
        opt_settings={"constraints": []},
    )

    workflow = stjames.ScanWorkflow(
        initial_molecule=initial_molecule,
        scan_settings=scan_settings,
        calc_settings=calc_settings,
        calc_engine=calculation_engine or calculation_method.default_engine(),
        wavefront_propagation=wavefront_propagation,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "scan",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": initial_molecule,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())
