"""Covalent inhibitor scan workflow - bond scan for a covalent inhibitor reaction."""

from dataclasses import dataclass

import stjames

from ..folder import Folder
from ..molecule import Molecule
from ..protein import Protein
from ..utils import api_client
from .base import Workflow, WorkflowResult, register_result


@dataclass(frozen=True, slots=True)
class CovalentInhibitorScanPoint:
    """
    A point along the covalent inhibitor bond scan.

    :param index: index of the point along the scan
    :param distance: scanned bond distance, in Å
    :param molecule: Molecule at the point
    :param energy: implicit-solvent single-point energy at this geometry, in Hartree
    :param uuid: UUID of the scan-point calculation
    """

    index: int
    distance: float
    molecule: Molecule
    energy: float | None
    uuid: str | None


@register_result("covalent_inhibitor_scan")
class CovalentInhibitorScanResult(WorkflowResult):
    """Result from a covalent inhibitor scan workflow."""

    _stjames_class = stjames.CovalentInhibitorScanWorkflow

    def __repr__(self) -> str:
        return f"<CovalentInhibitorScanResult points={len(self.scan_points)}>"

    @property
    def scan_points(self) -> list[CovalentInhibitorScanPoint]:
        """Points along the scan, ordered by scan index (scan_start to scan_stop)."""
        return [
            CovalentInhibitorScanPoint(
                index=p.index,
                distance=p.distance,
                molecule=Molecule.from_stjames(p.molecule),
                energy=p.energy,
                uuid=p.uuid,
            )
            for p in self._workflow.scan_points
        ]

    def get_energies(self) -> list[tuple[float, float | None]]:
        """
        Get scan distances paired with single-point energies.

        :returns: List of (distance, energy) tuples, in Å and Hartree respectively.
        """
        return [(p.distance, p.energy) for p in self.scan_points]


def submit_covalent_inhibitor_scan_workflow(
    protein: str | Protein,
    protein_reactive_atom_index: int,
    ligand_reactive_atom_index: int,
    settings: stjames.CovalentInhibitorScanSettings | None = None,
    name: str = "Covalent Inhibitor Scan Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
    is_draft: bool = False,
) -> Workflow:
    """
    Submits a covalent inhibitor scan workflow to the API.

    :param protein: covalently docked protein-ligand complex (protein plus the ligand as a
        non-polymer residue). Can be a UUID or a Protein object.
    :param protein_reactive_atom_index: 0-based index of the reacting protein atom, in PDB
        record order.
    :param ligand_reactive_atom_index: 0-based index of the reacting ligand atom, in PDB
        record order.
    :param settings: settings controlling the scan. Defaults to
        `stjames.CovalentInhibitorScanSettings()`.
    :param name: name of the workflow.
    :param folder_uuid: UUID of the folder to place the workflow in.
    :param folder: Folder object to store the workflow in.
    :param max_credits: maximum number of credits to use for the workflow.
    :param webhook_url: URL that Rowan will POST to when the workflow completes.
    :param is_draft: if True, submit the workflow as a draft without starting execution.
    :returns: Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder:
        folder_uuid = folder.uuid
    if isinstance(protein, Protein):
        protein = protein.uuid

    workflow = stjames.CovalentInhibitorScanWorkflow(
        protein=protein,
        protein_reactive_atom_index=protein_reactive_atom_index,
        ligand_reactive_atom_index=ligand_reactive_atom_index,
        settings=settings or stjames.CovalentInhibitorScanSettings(),
    )

    data = {
        "workflow_type": "covalent_inhibitor_scan",
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
