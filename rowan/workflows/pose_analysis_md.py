"""Pose-analysis MD workflow - molecular dynamics simulations for ligand-protein complexes."""

from dataclasses import dataclass
from pathlib import Path

import stjames

from ..protein import Protein, retrieve_protein
from ..utils import api_client
from .base import Message, Workflow, WorkflowResult, parse_messages, register_result


@dataclass(frozen=True, slots=True)
class TrajectoryResult:
    """Results from a single MD trajectory replicate."""

    uuid: str
    """UUID of the trajectory calculation."""

    ligand_rmsd: list[float]
    """Ligand RMSD values over time (Angstrom)."""

    contacts: dict
    """Contact analysis data between ligand and protein."""

    cluster_indices_by_frame: list[int]
    """Cluster assignment for each frame (0-indexed)."""

    cluster_centroid_indices: list[int]
    """Frame indices of cluster centroids."""


@register_result("pose_analysis_md")
class PoseAnalysisMDResult(WorkflowResult):
    """Result from a Pose-Analysis Molecular Dynamics (MD) workflow."""

    _stjames_class = stjames.PoseAnalysisMolecularDynamicsWorkflow

    def __repr__(self) -> str:
        n_traj = len(self.trajectories)
        n_clust = self.num_clusters
        return f"<PoseAnalysisMDResult trajectories={n_traj} clusters={n_clust}>"

    @property
    def num_clusters(self) -> int | None:
        """Number of conformational clusters identified across all trajectories."""
        return getattr(self._workflow, "num_clusters", None)

    @property
    def trajectories(self) -> list[TrajectoryResult]:
        """
        Results from each trajectory replicate.

        Each trajectory contains RMSD values, contact analysis, and cluster assignments.
        """
        raw = getattr(self._workflow, "trajectories", []) or []
        results: list[TrajectoryResult] = []
        for t in raw:
            results.append(
                TrajectoryResult(
                    uuid=t.uuid,
                    ligand_rmsd=list(t.ligand_rmsd or []),
                    contacts=dict(t.contacts) if t.contacts else {},
                    cluster_indices_by_frame=list(t.cluster_indices_by_frame or []),
                    cluster_centroid_indices=list(t.cluster_centroid_indices or []),
                )
            )
        return results

    @property
    def minimized_protein_uuid(self) -> str | None:
        """UUID of the energy-minimized protein structure."""
        return getattr(self._workflow, "minimized_protein_uuid", None)

    def get_minimized_protein(self) -> Protein | None:
        """
        Fetch the energy-minimized protein structure.

        :return: Protein object or None if not available.
        """
        uuid = self.minimized_protein_uuid
        if not uuid:
            return None
        if "minimized_protein" not in self._cache:
            self._cache["minimized_protein"] = retrieve_protein(uuid)
        return self._cache["minimized_protein"]

    @property
    def messages(self) -> list[Message]:
        """Any messages or warnings from the workflow."""
        return parse_messages(getattr(self._workflow, "messages", None))


    def download_trajectories(
        self,
        replicates: list[int],
        name: str | None = None,
        path: Path | None = None,
    ) -> Path:
        """
        Download DCD trajectory files for specified replicates.

        :param replicates: List of replicate indices to download.
        :param name: Custom name for the tar.gz file (without extension).
        :param path: Directory to save the file to. Defaults to current directory.
        :return: Path to the downloaded tar.gz file.
        :raises HTTPError: If the API request fails.
        """
        if path is None:
            path = Path.cwd()

        path.mkdir(parents=True, exist_ok=True)

        with api_client() as client:
            response = client.post(
                f"/trajectory/{self.workflow_uuid}/trajectory_dcds",
                json=replicates,
            )
            response.raise_for_status()

        file_name = f"{name or 'trajectories'}.tar.gz"
        file_path = path / file_name
        with open(file_path, "wb") as f:
            f.write(response.content)

        return file_path


def submit_pose_analysis_md_workflow(
    protein: str | Protein,
    initial_smiles: str,
    num_trajectories: int = 4,
    equilibration_time_ns: float = 1,
    simulation_time_ns: float = 10,
    temperature: float = 300,
    pressure_atm: float = 1.0,
    langevin_timescale_ps: float = 1.0,
    timestep_fs: float = 2,
    constrain_hydrogens: bool = True,
    nonbonded_cutoff: float = 8.0,
    ionic_strength_M: float = 0.10,
    water_buffer: float = 6.0,
    ligand_residue_name: str = "LIG",
    protein_restraint_cutoff: float = 7.0,
    protein_restraint_constant: float = 100,
    save_solvent: bool = False,
    name: str = "Pose-Analysis MD Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a Pose-Analysis Molecular Dynamics (MD) workflow to the API.

    :param protein: The *holo* protein on which MD will be run.
        Can be input as a UUID or a Protein object.
    :param initial_smiles: The SMILES for the ligand.
    :param num_trajectories: The number of trajectories to run.
    :param equilibration_time_ns: how long to equilibrate trajectories for, in ns
    :param simulation_time_ns: how long to run trajectories for, in ns
    :param temperature: temperature, in K
    :param pressure_atm: pressure, in atm
    :param langevin_timescale_ps: timescale for the Langevin integrator, in ps⁻¹
    :param timestep_fs: timestep, in femtoseconds
    :param ligand_residue_name: The name of the residue corresponding to the ligand.
    :param constrain_hydrogens: whether or not to use SHAKE to freeze bonds to hydrogen
    :param nonbonded_cutoff: nonbonded cutoff for particle-mesh Ewald, in Å
    :param ionic_strength_M: ionic strength of the solution, in M (molar)
    :param water_buffer: amount of water to add around the protein, in Å
    :param protein_restraint_cutoff: cutoff past which alpha-carbons will be constrained, in Å
    :param protein_restraint_constant: force constant for backbone restraints, in kcal/mol/Å²
    :param save_solvent: whether to save solvent molecules
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if isinstance(protein, Protein):
        protein = protein.uuid

    workflow = stjames.PoseAnalysisMolecularDynamicsWorkflow(
        protein=protein,
        initial_smiles=initial_smiles,
        num_trajectories=num_trajectories,
        equilibration_time_ns=equilibration_time_ns,
        simulation_time_ns=simulation_time_ns,
        temperature=temperature,
        pressure_atm=pressure_atm,
        langevin_timescale_ps=langevin_timescale_ps,
        timestep_fs=timestep_fs,
        ligand_residue_name=ligand_residue_name,
        constrain_hydrogens=constrain_hydrogens,
        nonbonded_cutoff=nonbonded_cutoff,
        ionic_strength_M=ionic_strength_M,
        water_buffer=water_buffer,
        protein_restraint_cutoff=protein_restraint_cutoff,
        protein_restraint_constant=protein_restraint_constant,
        save_solvent=save_solvent,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "pose_analysis_md",
        "workflow_data": workflow.model_dump(serialize_as_any=True, mode="json"),
        "initial_smiles": initial_smiles,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


__all__ = ["PoseAnalysisMDResult", "TrajectoryResult", "submit_pose_analysis_md_workflow"]
