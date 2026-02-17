"""Scan workflow - perform potential energy surface scans."""

from typing import Any

import stjames
from rdkit import Chem

from ..utils import api_client
from .base import RdkitMol, StJamesMolecule, Workflow, WorkflowResult, register_result


@register_result("scan")
class ScanResult(WorkflowResult):
    """Result from a scan workflow."""

    _stjames_class = stjames.ScanWorkflow

    def __repr__(self) -> str:
        n = len(self.scan_points)
        return f"<ScanResult scan_points={n}>"

    @property
    def scan_points(self) -> list[str | None]:
        """UUIDs of scan point calculations."""
        return list(self._workflow.scan_points)


def submit_scan_workflow(
    initial_molecule: dict[str, Any] | StJamesMolecule | RdkitMol,
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

    :param initial_molecule: The molecule used in the scan.
    :param scan_settings: The scan settings.
    :param calculation_engine: The engine to use for the calculation.
    :param calculation_method: The method to use for the calculation.
    :param wavefront_propagation: Whether to use wavefront propagation in the scan.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to store the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if isinstance(initial_molecule, StJamesMolecule):
        initial_molecule = initial_molecule.model_dump(mode="json")
    elif isinstance(initial_molecule, Chem.rdchem.Mol | Chem.rdchem.RWMol):
        initial_molecule = StJamesMolecule.from_rdkit(initial_molecule, cid=0)

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


__all__ = ["ScanResult", "submit_scan_workflow"]
