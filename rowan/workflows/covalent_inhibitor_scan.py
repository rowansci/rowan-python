"""Covalent inhibitor scan workflow - reactivity profile via ML/MM umbrella sampling."""

from dataclasses import dataclass

import stjames

from ..folder import Folder
from ..molecule import Molecule
from ..protein import Protein
from ..types import ProteinUUID
from ..utils import api_client
from .base import Workflow, WorkflowResult, register_result


@dataclass(frozen=True, slots=True)
class CovalentInhibitorScanPoint:
    """
    A point on the covalent inhibitor reactivity profile.

    :param index: window index.
    :param distance: bias center for this window, in Å.
    :param force_constant: harmonic bias force constant, in kcal/mol/Å².
    :param free_energy: unbiased free energy at this center, in kcal/mol, if available.
    :param mean_distance: mean sampled reactive-bond distance, in Å, if available.
    :param n_samples: production samples contributing to this window.
    :param molecule: final sampled geometry of the model region, if available.
    """

    index: int
    distance: float
    force_constant: float
    free_energy: float | None
    mean_distance: float | None
    n_samples: int
    molecule: Molecule | None


@dataclass(frozen=True, slots=True)
class UmbrellaSamplingConvergence:
    """
    Diagnostics saying whether a reactivity profile can be trusted.

    :param overlap_matrix: MBAR overlap between windows; low overlap between neighbours
        means the profile should not be trusted.
    :param round_trips: replicas that traversed the whole coordinate and returned.
    :param worst_pair_acceptance: lowest exchange acceptance over neighbouring pairs, if available.
    """

    overlap_matrix: tuple[tuple[float, ...], ...]
    round_trips: int
    worst_pair_acceptance: float | None


@register_result("covalent_inhibitor_scan")
class CovalentInhibitorScanResult(WorkflowResult):
    """Result from a covalent inhibitor scan workflow."""

    _stjames_class = stjames.CovalentInhibitorScanWorkflow

    def __repr__(self) -> str:
        return f"<CovalentInhibitorScanResult points={len(self.points)} barrier={self.barrier}>"

    @property
    def points(self) -> list[CovalentInhibitorScanPoint]:
        """Points on the reactivity profile, in window-index order."""
        return [
            CovalentInhibitorScanPoint(
                index=p.index,
                distance=p.distance,
                force_constant=p.force_constant,
                free_energy=p.free_energy,
                mean_distance=p.mean_distance,
                n_samples=p.n_samples,
                molecule=Molecule.from_stjames(p.molecule) if p.molecule is not None else None,
            )
            for p in self._workflow.points
        ]

    @property
    def convergence(self) -> UmbrellaSamplingConvergence | None:
        """Sampling convergence diagnostics, when available."""
        c = self._workflow.convergence
        if c is None:
            return None
        return UmbrellaSamplingConvergence(
            overlap_matrix=tuple(tuple(row) for row in c.overlap_matrix),
            round_trips=c.round_trips,
            worst_pair_acceptance=c.worst_pair_acceptance,
        )

    @property
    def barrier(self) -> float | None:
        """Activation free energy for addition, in kcal/mol, or None without an interior maximum."""
        return self._workflow.barrier

    def get_energies(self) -> list[tuple[float, float | None]]:
        """
        Get bias-center distances paired with free energies.

        :returns: List of (distance, free_energy) tuples, in Å and kcal/mol respectively.
        """
        return [(p.distance, p.free_energy) for p in self.points]


def submit_covalent_inhibitor_scan_workflow(
    protein: Protein | ProteinUUID,
    protein_reactive_atom_index: int,
    ligand_reactive_atom_index: int,
    reactant_smiles: str,
    settings: stjames.UmbrellaSamplingScanSettings | None = None,
    name: str = "Covalent Inhibitor Scan Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
    is_draft: bool = False,
) -> Workflow:
    """
    Submits a covalent inhibitor scan workflow to the API.

    .. warning::
        This workflow is in beta. Its interface and behavior may change.

    Seeds windows along the reactive bond distance with a steered pull, samples them via
    ML/MM umbrella sampling, and combines them into a free energy profile.

    :param protein: covalently docked protein-ligand complex (protein plus the ligand as a
        non-polymer residue). Can be a UUID or a Protein object.
    :param protein_reactive_atom_index: 0-based index of the reacting protein atom, in PDB
        record order.
    :param ligand_reactive_atom_index: 0-based index of the reacting ligand atom, in PDB
        record order.
    :param reactant_smiles: SMILES of the neutral reactant ligand, used to rebuild the ligand
        as a separate non-covalent molecule at the MM level.
    :param settings: settings controlling the umbrella sampling. Defaults to
        `stjames.UmbrellaSamplingScanSettings()`.
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
        reactant_smiles=reactant_smiles,
        settings=settings or stjames.UmbrellaSamplingScanSettings(),
    )

    data = {
        "workflow_type": "covalent_inhibitor_scan",
        "workflow_data": workflow.model_dump(serialize_as_any=True, mode="json"),
        "initial_smiles": reactant_smiles,
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
