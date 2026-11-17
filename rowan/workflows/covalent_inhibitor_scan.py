"""Covalent reaction free-energy surfaces from ML/MM umbrella sampling."""

from pathlib import Path

from rowan.folder import Folder
from rowan.protein import Protein, retrieve_protein
from rowan.types import ProteinUUID
from rowan.utils import api_client, download_presigned_files
from rowan.workflows.base import Message, Workflow, WorkflowResult, parse_messages, register_result
from stjames.workflows.covalent_inhibitor_scan import (
    CovalentInhibitorScanSettings,
    CovalentInhibitorScanWorkflow,
    FreeEnergySurface,
    ResolvedScanSetup,
    SurfaceDiagnostics,
    UmbrellaWindowResult,
)


@register_result("covalent_inhibitor_scan")
class CovalentInhibitorScanResult(WorkflowResult):
    """Reaction-coordinate windows, free-energy surface, and sampling diagnostics."""

    _stjames_class = CovalentInhibitorScanWorkflow

    def __repr__(self) -> str:
        return (
            f"<CovalentInhibitorScanResult windows={len(self.windows)} "
            f"surface={self.surface is not None}>"
        )

    @property
    def settings(self) -> CovalentInhibitorScanSettings:
        """Submitted reactive region, scan coordinates, and simulation settings."""
        return self._workflow.settings

    @property
    def resolved_setup(self) -> ResolvedScanSetup | None:
        """Input atom mappings, capped bonds, and reactive-region charge and spin."""
        return getattr(self._workflow, "resolved_setup", None)

    @property
    def windows(self) -> list[UmbrellaWindowResult]:
        """Window results in submitted grid order, including snapshots and DCD references.

        Each `molecule` is a `rowan.StJamesMolecule`, or `None`. Use
        `rowan.Molecule.from_stjames` for SDK molecule helpers.
        """
        return getattr(self._workflow, "windows", [])

    @property
    def surface(self) -> FreeEnergySurface | None:
        """Free energies and uncertainties, in kcal/mol.

        Flattened arrays follow ascending bin edges, with the last axis varying fastest.
        Missing cell values remain `None`.
        """
        return getattr(self._workflow, "surface", None)

    @property
    def surface_diagnostics(self) -> SurfaceDiagnostics | None:
        """Exchange counts and MBAR overlap in window-list order."""
        return getattr(self._workflow, "surface_diagnostics", None)

    @property
    def topology_uuid(self) -> str | None:
        """UUID of the shared solvated topology in DCD atom order."""
        return getattr(self._workflow, "topology_uuid", None)

    @property
    def messages(self) -> list[Message]:
        """Workflow warnings, errors, and informational messages."""
        return parse_messages(getattr(self._workflow, "messages", None))

    def get_topology(self) -> Protein | None:
        """Fetch and cache the full-system topology shared by the window trajectories.

        Returns:
            Topology protein, or `None` before a topology is saved
        """
        if not (uuid := self.topology_uuid):
            return None
        if "topology" not in self._cache:
            self._cache["topology"] = retrieve_protein(uuid, workflow_uuid=self.workflow_uuid)
        return self._cache["topology"]

    def download_trajectory_files(
        self,
        window_indices: list[int],
        path: Path | str | None = None,
    ) -> list[Path]:
        """Download full-system DCDs for selected positions in `windows`.

        Args:
            window_indices: zero-based positions in the window list, not grid index tuples
            path: output directory, defaulting to the current directory

        Returns:
            Downloaded paths in requested window order

        Raises:
            IndexError: window position is negative or outside the window list
            ValueError: selected window has no saved trajectory
            httpx.HTTPStatusError: API or storage request fails
        """
        windows = self.windows
        for index in window_indices:
            if index < 0 or index >= len(windows):
                raise IndexError(f"Window position {index} out of range (0-{len(windows) - 1})")
            if windows[index].trajectory is None:
                raise ValueError(f"Window at position {index} has no saved trajectory")
        if not window_indices:
            return []
        return download_presigned_files(
            Path(path) if path is not None else Path.cwd(),
            f"/trajectory/{self.workflow_uuid}/trajectory_dcd_urls",
            json=window_indices,
        )


def submit_covalent_inhibitor_scan_workflow(
    protein: Protein | ProteinUUID,
    reactant_smiles: str,
    settings: CovalentInhibitorScanSettings,
    name: str = "Covalent Inhibitor Scan Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
    is_draft: bool = False,
) -> Workflow[CovalentInhibitorScanResult]:
    """Submit a covalent reaction free-energy scan in a prepared protein–ligand complex.

    Warning:
        This workflow is in beta. Its interface and behavior may change.

    One or two distance coordinates, including weighted distance combinations, define
    the umbrella grid. Atom references use residue identifiers and input atom names.
    The caller supplies explicit hydrogens and the intended protonation and connectivity.

    Args:
        protein: prepared protein–ligand complex, as a UUID or Protein object
        reactant_smiles: ligand reactant SMILES matching the intended protonation and connectivity
        settings: ligand selection, reactive region, scan coordinates, and simulation settings
        name: name of the workflow
        folder_uuid: UUID of the folder to place the workflow in
        folder: destination folder
        max_credits: maximum credits for the workflow
        webhook_url: URL that Rowan will POST to when the workflow completes
        is_draft: save as a draft without starting execution

    Returns:
        Submitted workflow with typed surface results

    Raises:
        ValueError: both folder arguments are supplied, or the shared schema rejects the inputs
        httpx.HTTPStatusError: request to the API fails
    """
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder:
        folder_uuid = folder.uuid
    if isinstance(protein, Protein):
        protein = protein.uuid

    workflow = CovalentInhibitorScanWorkflow(
        protein=protein,
        reactant_smiles=reactant_smiles,
        settings=settings,
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
        return Workflow[CovalentInhibitorScanResult](**response.json())
