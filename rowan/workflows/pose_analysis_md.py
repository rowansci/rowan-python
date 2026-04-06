"""Pose-analysis MD workflow - molecular dynamics simulations for ligand-protein complexes."""

from dataclasses import dataclass
from pathlib import Path

import stjames

from ..folder import Folder
from ..protein import Protein, retrieve_protein
from ..utils import api_client
from .base import Message, Workflow, WorkflowResult, parse_messages, register_result


@dataclass(frozen=True, slots=True)
class TrajectoryResult:
    """
    Results from a single MD trajectory replicate.

    :param uuid: UUID of the trajectory calculation.
    :param ligand_rmsd: Ligand RMSD values over time (Angstrom).
    :param contacts: Ligand-protein contacts with occupancy over the trajectory.
    """

    uuid: str
    ligand_rmsd: list[float]
    contacts: list[stjames.BindingPoseContact]


@register_result("pose_analysis_md")
class PoseAnalysisMDResult(WorkflowResult):
    """Result from a Pose-Analysis Molecular Dynamics (MD) workflow."""

    _stjames_class = stjames.PoseAnalysisMolecularDynamicsWorkflow

    def __repr__(self) -> str:
        n_traj = len(self.trajectories)
        return f"<PoseAnalysisMDResult trajectories={n_traj}>"

    @property
    def trajectories(self) -> list[TrajectoryResult]:
        """
        Results from each trajectory replicate.

        Each trajectory contains RMSD values, contact analysis, and cluster assignments.
        """
        raw = self._workflow.trajectories or []
        return [
            TrajectoryResult(uuid=t.uuid, ligand_rmsd=t.ligand_rmsd, contacts=t.contacts)
            for t in raw
        ]

    @property
    def average_rmsds(self) -> list[float | None]:
        """Average ligand RMSD per trajectory (Angstrom)."""
        return [
            sum(t.ligand_rmsd) / len(t.ligand_rmsd) if t.ligand_rmsd else None
            for t in self.trajectories
        ]

    @property
    def minimized_protein_uuid(self) -> str | None:
        """UUID of the energy-minimized protein structure."""
        return getattr(self._workflow, "minimized_protein_uuid", None)

    def get_minimized_protein(self) -> Protein | None:
        """
        Fetch the energy-minimized protein structure.

        .. note::
            Makes one API call on first access.
            Results are cached. Call clear_cache() to refresh.

        :returns: Protein object or None if not available.
        """
        if not (uuid := self.minimized_protein_uuid):
            return None
        if "minimized_protein" not in self._cache:
            self._cache["minimized_protein"] = retrieve_protein(uuid)
        return self._cache["minimized_protein"]

    @property
    def messages(self) -> list[Message]:
        """Any messages or warnings from the workflow."""
        return parse_messages(self._workflow.messages)

    def get_atom_distances(
        self,
        atom_pairs: list[tuple[int, int]],
        replicate: int = 0,
    ) -> list[list[float]]:
        """
        Fetch interatomic distances over the trajectory for specified atom pairs.

        Atom indices can be found in the ``contacts`` field of each trajectory,
        which provides ``ligand_atom_index`` and ``protein_atom_index`` for each contact.

        :param atom_pairs: List of (atom_i, atom_j) index pairs (0-indexed).
        :param replicate: Trajectory replicate index (default 0).
        :returns: List of distance arrays, one per pair, over all frames (Angstrom).
        :raises HTTPError: If the API request fails.
        """
        with api_client() as client:
            response = client.post(
                f"/trajectory/{self.workflow_uuid}/atom_trajectories",
                params={"replicate": replicate},
                json=atom_pairs,
            )
            response.raise_for_status()
        return response.json()

    def download_trajectories(
        self,
        replicates: list[int],
        name: str | None = None,
        path: Path | str | None = None,
    ) -> Path:
        """
        Download DCD trajectory files for specified replicates.

        :param replicates: List of replicate indices to download.
        :param name: Custom name for the tar.gz file (without extension).
        :param path: Directory to save the file to. Defaults to current directory.
        :returns: Path to the downloaded tar.gz file.
        :raises HTTPError: If the API request fails.
        """
        path = Path(path) if path is not None else Path.cwd()
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
    validate_forcefield: bool = True,
    name: str = "Pose-Analysis MD Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
    is_draft: bool = False,
) -> Workflow:
    """
    Submits a Pose-Analysis Molecular Dynamics (MD) workflow to the API.

    :param protein: *Holo* protein on which MD will be run.
        Can be input as a UUID or a Protein object.
    :param initial_smiles: SMILES for the ligand.
    :param num_trajectories: Number of trajectories to run.
    :param equilibration_time_ns: Equilibration time per trajectory, in ns.
    :param simulation_time_ns: Simulation time per trajectory, in ns.
    :param temperature: Temperature, in K.
    :param pressure_atm: Pressure, in atm.
    :param langevin_timescale_ps: Timescale for the Langevin integrator, in ps⁻¹.
    :param timestep_fs: Timestep, in femtoseconds.
    :param ligand_residue_name: Name of the residue corresponding to the ligand.
    :param constrain_hydrogens: Whether to use SHAKE to freeze bonds to hydrogen.
    :param nonbonded_cutoff: Nonbonded cutoff for particle-mesh Ewald, in Å.
    :param ionic_strength_M: Ionic strength of the solution, in M (molar).
    :param water_buffer: Amount of water to add around the protein, in Å.
    :param protein_restraint_cutoff: Cutoff past which alpha-carbons will be constrained, in Å.
    :param protein_restraint_constant: Force constant for backbone restraints, in kcal/mol/Å².
    :param save_solvent: Whether to save solvent molecules.
    :param validate_forcefield: if True (default), validate the protein forcefield
        compatibility before submitting. Raises an error early if the protein cannot
        be parameterized or has clashing residues.
    :param name: Name of the workflow.
    :param folder_uuid: UUID of the folder to place the workflow in.
    :param folder: Folder object to store the workflow in.
    :param max_credits: Maximum number of credits to use for the workflow.
    :param webhook_url: URL that Rowan will POST to when the workflow completes.
    :param is_draft: If True, submit the workflow as a draft without starting execution.
    :returns: Workflow object representing the submitted workflow.
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
        "workflow_type": "pose_analysis_md",
        "workflow_data": workflow.model_dump(serialize_as_any=True, mode="json"),
        "initial_smiles": initial_smiles,
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
