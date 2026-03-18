"""RBFE perturbation workflow - run relative binding free energy FEP simulations."""

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import stjames
from stjames.workflows.relative_binding_free_energy_perturbation import TMDRBFESettings

from ..folder import Folder
from ..molecule import Molecule
from ..protein import Protein
from ..utils import api_client
from .base import Workflow, WorkflowResult, molecule_to_dict, register_result
from .rbfe_graph import RelativeBindingFreeEnergyGraphEdge, RelativeBindingFreeEnergyGraphResult

_PRESETS: dict[str, dict] = {
    "fast": {
        "n_eq_steps": 10_000,
        "n_frames": 1_000,
        "n_windows": 24,
        "min_overlap": 0.2,
        "target_overlap": 0.2,
        "local_md_steps": 390,
        "charge_method": "nagl",
    },
    "recommended": {
        "n_eq_steps": 200_000,
        "n_frames": 2_000,
        "n_windows": 48,
        "min_overlap": 0.667,
        "target_overlap": 0.667,
        "local_md_steps": 390,
        "charge_method": "nagl",
    },
    "rigorous": {
        "n_eq_steps": 200_000,
        "n_frames": 2_000,
        "n_windows": 48,
        "min_overlap": 0.667,
        "target_overlap": 0.667,
        "local_md_steps": 0,
        "charge_method": "amber_am1bcc",
    },
}


@dataclass(frozen=True, slots=True)
class RelativeBindingFreeEnergyResult:
    """Aggregate RBFE outcome for a single ligand.

    :param dg: Predicted binding free energy difference (kcal/mol).
    :param dg_err: Uncertainty estimate on dg.
    """

    dg: float
    dg_err: float


@dataclass(frozen=True, slots=True)
class RelativeBindingFreeEnergyDiagnostics:
    """Quality-control metrics from an RBFE simulation.

    :param cycle_closure_rms: RMS error across completed thermodynamic cycles.
    :param windows_completed: Count of successfully converged lambda windows.
    :param windows_failed: Count of failed lambda windows.
    """

    cycle_closure_rms: float | None
    windows_completed: int | None
    windows_failed: int | None


@register_result("relative_binding_free_energy_perturbation")
class RelativeBindingFreeEnergyPerturbationResult(WorkflowResult):
    """Result from a relative binding free energy perturbation workflow."""

    _stjames_class = stjames.RelativeBindingFreeEnergyPerturbationWorkflow

    def __repr__(self) -> str:
        n = len(self.ligand_dg_results) if self.ligand_dg_results else 0
        return f"<RelativeBindingFreeEnergyPerturbationResult ligands={n}>"

    @property
    def ligands(self) -> dict[str, Molecule]:
        """Ligand molecules keyed by identifier."""
        return {k: Molecule.from_stjames(v) for k, v in self._workflow.ligands.items()}

    @property
    def edges(self) -> list[RelativeBindingFreeEnergyGraphEdge]:
        """Graph edges with per-edge FEP results."""
        g = self._workflow.graph
        return [
            RelativeBindingFreeEnergyGraphEdge(
                ligand_a=e.ligand_a,
                ligand_b=e.ligand_b,
                core=e.core,
                score=e.score,
                ddg=e.ddg,
                ddg_err=e.ddg_err,
                complex_dg=e.complex_dg,
                complex_dg_err=e.complex_dg_err,
                solvent_dg=e.solvent_dg,
                solvent_dg_err=e.solvent_dg_err,
                vacuum_dg=e.vacuum_dg,
                vacuum_dg_err=e.vacuum_dg_err,
                failed=e.failed,
                complex_lambda_values=e.complex_lambda_values,
                complex_overlap_matrix=e.complex_overlap_matrix,
            )
            for e in g.edges
        ]

    def download_edge_trajectories(
        self,
        edge_index: int,
        path: Path | str | None = None,
        name: str | None = None,
    ) -> Path:
        """
        Download all DCD trajectory files for a specific perturbation edge.

        :param edge_index: Index of the edge (0-based, matching ``edges`` order).
        :param path: Directory to save the file to. Defaults to current directory.
        :param name: Custom name for the tar.gz file (without extension).
        :returns: Path to the downloaded tar.gz file.
        :raises IndexError: If edge_index is out of range.
        :raises HTTPError: If the API request fails.
        """
        edges = self.edges
        if edge_index < 0 or edge_index >= len(edges):
            raise IndexError(f"Edge index {edge_index} out of range (0-{len(edges) - 1})")

        path = Path(path) if path is not None else Path.cwd()
        path.mkdir(parents=True, exist_ok=True)

        with api_client() as client:
            response = client.post(
                f"/trajectory/{self.workflow_uuid}/rbfe_trajectory_dcds",
                params={"edge_index": edge_index},
            )
            response.raise_for_status()

        file_name = f"{name or f'edge_{edge_index}_trajectories'}.tar.gz"
        file_path = path / file_name
        with open(file_path, "wb") as f:
            f.write(response.content)

        return file_path

    def download_all_trajectories(
        self,
        path: Path | str | None = None,
    ) -> list[Path]:
        """
        Download DCD trajectory files for all perturbation edges.

        :param path: Directory to save the files to. Defaults to current directory.
        :returns: List of paths to the downloaded tar.gz files, one per edge.
        :raises HTTPError: If any API request fails.
        """
        return [self.download_edge_trajectories(i, path=path) for i in range(len(self.edges))]

    @property
    def ligand_dg_results(self) -> dict[str, RelativeBindingFreeEnergyResult] | None:
        """Per-ligand binding free energy results, or None if not yet computed."""
        raw = self._workflow.ligand_dg_results
        if raw is None:
            return None
        return {
            k: RelativeBindingFreeEnergyResult(dg=v.dg, dg_err=v.dg_err) for k, v in raw.items()
        }

    @property
    def diagnostics(self) -> RelativeBindingFreeEnergyDiagnostics | None:
        """Aggregate QC metrics from the FEP simulation."""
        d = self._workflow.diagnostics
        if d is None:
            return None
        return RelativeBindingFreeEnergyDiagnostics(
            cycle_closure_rms=d.cycle_closure_rms,
            windows_completed=d.windows_completed,
            windows_failed=d.windows_failed,
        )


def submit_relative_binding_free_energy_perturbation_workflow(
    graph_result: RelativeBindingFreeEnergyGraphResult,
    protein: str | Protein,
    tmd_settings: Literal["fast", "recommended", "rigorous"] = "recommended",
    forcefield: Literal["off_sage_2_0_0", "off_sage_2_2_1"] = "off_sage_2_0_0",
    charge_method: Literal["amber_am1bcc", "nagl"] | None = None,
    n_eq_steps: int | None = None,
    n_frames: int | None = None,
    steps_per_frame: int = 400,
    n_windows: int | None = None,
    min_overlap: float | None = None,
    target_overlap: float | None = None,
    water_sampling_padding: float = 0.4,
    rest_max_temperature_scale: float = 1.0,
    rest_temperature_scale_interpolation: Literal["exponential", "linear"] = "exponential",
    local_md_steps: int | None = None,
    local_md_k: float = 10_000.0,
    local_md_radius: float = 1.2,
    local_md_free_reference: bool = False,
    legs: list[Literal["vacuum", "solvent", "complex"]] | None = None,
    save_trajectories: bool = False,
    trajectory_save_interval: int = 1000,
    validate_forcefield: bool = True,
    name: str = "RBFE Perturbation",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a relative binding free energy perturbation (RBFE) workflow to the API.

    Runs FEP simulations along edges of the perturbation graph to predict
    relative binding free energies between ligands.

    Preset settings (any individual param overrides the tmd_settings):
    - ``"fast"``: fewer windows/steps for quick screening (NAGL charges).
    - ``"recommended"`` (default): balanced speed and accuracy (NAGL charges).
    - ``"rigorous"``: same as recommended but disables local MD for higher accuracy.

    :param graph_result: Completed ``RelativeBindingFreeEnergyGraphResult``.
    :param protein: Protein target, as a UUID string or Protein object.
    :param tmd_settings: Starting settings profile. Individual params override this.
    :param forcefield: Force field for the simulation (e.g. ``"off_sage_2_0_0"``).
    :param charge_method: Method for computing partial charges.
    :param n_eq_steps: Equilibration steps per lambda window.
    :param n_frames: Production frames saved per lambda window.
    :param steps_per_frame: MD integration steps per saved frame.
    :param n_windows: Maximum number of lambda windows considered for bisection.
    :param min_overlap: Minimum acceptable overlap during schedule bisection.
    :param target_overlap: Desired overlap after HREX optimization.
    :param water_sampling_padding: Extra nanometers added to the solvent sampling radius.
    :param rest_max_temperature_scale: Maximum effective temperature scaling for REST.
    :param rest_temperature_scale_interpolation: Functional form used for REST scaling.
    :param local_md_steps: Number of local MD steps per frame (0 disables local MD).
    :param local_md_k: Spring constant used during local MD.
    :param local_md_radius: Sphere radius in nanometers for the local MD region.
    :param local_md_free_reference: Whether to free the reference frame during local MD.
    :param legs: Which thermodynamic cycle legs to run.
    :param save_trajectories: Whether to save DCD trajectories.
    :param trajectory_save_interval: Save every Nth frame when saving trajectories.
    :param validate_forcefield: If True (default), validate protein forcefield
        compatibility before submitting.
    :param name: Name of the workflow.
    :param folder_uuid: UUID of the folder to place the workflow in.
    :param folder: Folder object to store the workflow in.
    :param max_credits: Maximum number of credits to use for the workflow.
    :returns: Workflow object representing the submitted workflow.
    :raises ValueError: If graph_result has no graph or both folder and folder_uuid are provided.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder:
        folder_uuid = folder.uuid

    if isinstance(protein, Protein):
        protein = protein.uuid

    if validate_forcefield:
        Protein(uuid=protein).validate_protein_forcefield()

    if graph_result.graph is None:
        raise ValueError("RBFEGraphResult has no graph — has the workflow completed successfully?")

    # Apply tmd_settings, then override with any explicitly provided params
    p = _PRESETS[tmd_settings]
    settings = TMDRBFESettings(
        forcefield=forcefield,
        charge_method=charge_method if charge_method is not None else p["charge_method"],
        n_eq_steps=n_eq_steps if n_eq_steps is not None else p["n_eq_steps"],
        n_frames=n_frames if n_frames is not None else p["n_frames"],
        steps_per_frame=steps_per_frame,
        n_windows=n_windows if n_windows is not None else p["n_windows"],
        min_overlap=min_overlap if min_overlap is not None else p["min_overlap"],
        target_overlap=target_overlap if target_overlap is not None else p["target_overlap"],
        water_sampling_padding=water_sampling_padding,
        rest_max_temperature_scale=rest_max_temperature_scale,
        rest_temperature_scale_interpolation=rest_temperature_scale_interpolation,
        local_md_steps=local_md_steps if local_md_steps is not None else p["local_md_steps"],
        local_md_k=local_md_k,
        local_md_radius=local_md_radius,
        local_md_free_reference=local_md_free_reference,
        legs=legs if legs is not None else ["solvent", "complex"],
        save_trajectories=save_trajectories,
        trajectory_save_interval=trajectory_save_interval,
    )

    ligands_dict = {k: molecule_to_dict(v) for k, v in graph_result.ligands.items()}

    workflow = stjames.RelativeBindingFreeEnergyPerturbationWorkflow(
        ligands=ligands_dict,
        graph=graph_result.graph,
        protein=protein,
        target=protein,
        settings=settings,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "relative_binding_free_energy_perturbation",
        "workflow_data": workflow.model_dump(mode="json"),
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())
