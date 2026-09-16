"""Pose-analysis MD workflow - molecular dynamics simulations for ligand-protein complexes."""

from dataclasses import dataclass, field
from typing import Literal

import stjames
from stjames import (
    GreedyClusteringSettings,
    KMeansClusteringSettings,
    ProteinForceField,
    WaterForceField,
)

from rowan.folder import Folder
from rowan.protein import Protein
from rowan.types import ProteinUUID
from rowan.utils import api_client

from ._molecular_dynamics import _MolecularDynamicsResult
from .base import Message, Workflow, parse_messages, register_result


@dataclass(frozen=True, slots=True)
class TrajectoryResult:
    """Results from a single MD trajectory replicate.

    Attributes:
        uuid: UUID of the trajectory calculation
        ligand_rmsd: ligand RMSD values over time (Angstrom)
        contacts: ligand-protein contacts with occupancy over the trajectory
        sasa: solvent-accessible surface area per analyzed frame (populated when
            analysis_interval_ps is set)
        polar_sasa: polar solvent-accessible surface area per analyzed frame (populated when
            analysis_interval_ps is set)
        isotropic_radius_of_gyration: radius of gyration per analyzed frame
        cluster_centroid_indices: frame indices of the cluster centroids (populated when
            clustering is set)
        cluster_indices_by_frame: cluster assignment for each frame (populated when clustering
            is set)
        protein_rmsd: per-frame Cα RMSD from the first frame, in angstrom
        rmsf: per-Cα RMSF from the mean structure, in angstrom
        potential_energy: per-frame potential energy of the simulated system, in Hartree
        mmgbsa_scores: per-frame MM/GBSA interaction energy, in kcal/mol
        mean_structure_uuid: UUID of the coordinate-averaged structure
        median_structure_frame_index: frame index of the medoid structure
    """

    uuid: str
    ligand_rmsd: list[float]
    contacts: list[stjames.BindingPoseContact]
    sasa: list[float | None]
    polar_sasa: list[float | None]
    isotropic_radius_of_gyration: list[float]
    cluster_centroid_indices: list[int]
    cluster_indices_by_frame: list[int]
    protein_rmsd: list[float] = field(default_factory=list)
    rmsf: list[float] = field(default_factory=list)
    potential_energy: list[float] = field(default_factory=list)
    mmgbsa_scores: list[float | None] = field(default_factory=list)
    mean_structure_uuid: str | None = None
    median_structure_frame_index: int | None = None


@register_result("pose_analysis_md")
class PoseAnalysisMDResult(_MolecularDynamicsResult):
    """Result from a Pose-Analysis Molecular Dynamics (MD) workflow."""

    _stjames_class = stjames.PoseAnalysisMolecularDynamicsWorkflow

    def __repr__(self) -> str:
        n_traj = len(self.trajectories)
        return f"<PoseAnalysisMDResult trajectories={n_traj}>"

    @property
    def trajectories(self) -> list[TrajectoryResult]:
        """Results from each trajectory replicate.

        Each trajectory contains RMSD values, contact analysis, and cluster assignments.
        """
        raw = self._workflow.trajectories or []
        return [
            TrajectoryResult(
                uuid=t.uuid,
                ligand_rmsd=t.binder_rmsd,
                contacts=t.contacts,
                sasa=t.sasa,
                polar_sasa=t.polar_sasa,
                isotropic_radius_of_gyration=t.isotropic_radius_of_gyration,
                cluster_centroid_indices=t.cluster_centroid_indices,
                cluster_indices_by_frame=t.cluster_indices_by_frame,
                protein_rmsd=t.protein_rmsd,
                rmsf=t.rmsf,
                potential_energy=t.potential_energy,
                mmgbsa_scores=t.mmgbsa_scores,
                mean_structure_uuid=t.mean_structure_uuid,
                median_structure_frame_index=t.median_structure_frame_index,
            )
            for t in raw
        ]

    @property
    def hydration_sites(self) -> list[stjames.HydrationSite]:
        """Hydration sites identified across all trajectories."""
        return self._workflow.hydration_sites or []

    @property
    def average_rmsds(self) -> list[float | None]:
        """Average ligand RMSD per trajectory (Angstrom)."""
        return [
            sum(t.ligand_rmsd) / len(t.ligand_rmsd) if t.ligand_rmsd else None
            for t in self.trajectories
        ]

    @property
    def messages(self) -> list[Message]:
        """Any messages or warnings from the workflow."""
        return parse_messages(self._workflow.messages)


def submit_pose_analysis_md_workflow(
    protein: Protein | ProteinUUID,
    initial_smiles: str,
    num_trajectories: int = 4,
    small_molecule_ff: Literal[
        "off_sage_2_0_0", "off_sage_2_2_1", "off_sage_2_3_0", "mango_1_0_0"
    ] = "off_sage_2_3_0",
    protein_ff: ProteinForceField | str = ProteinForceField.FF14SB,
    water_ff: WaterForceField | str = WaterForceField.TIP3P,
    equilibration_time_ns: float = 0.5,
    simulation_time_ns: float = 10,
    temperature: float = 300,
    pressure_atm: float = 1.0,
    langevin_timescale_ps: float = 1.0,
    timestep_fs: float = 4,
    hydrogen_mass: float = 3,
    constrain_hydrogens: bool = True,
    nonbonded_cutoff: float = 8.0,
    ionic_strength_M: float = 0.0,
    water_buffer: float = 8.0,
    ligand_residue_name: str = "LIG",
    protein_restraint_cutoff: float | None = 7.0,
    protein_restraint_constant: float = 100,
    save_solvent: bool = False,
    num_solvent_to_save: int | None = None,
    analysis_interval_ps: float | None = None,
    clustering: KMeansClusteringSettings | GreedyClusteringSettings | None = None,
    validate_forcefield: bool = True,
    name: str = "Pose-Analysis MD Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
    is_draft: bool = False,
) -> Workflow[PoseAnalysisMDResult]:
    """Submits a Pose-Analysis Molecular Dynamics (MD) workflow to the API.

    Args:
        protein: *Holo* protein on which MD will be run.
            Can be input as a UUID or a Protein object
        initial_smiles: SMILES for the ligand
        num_trajectories: number of trajectories to run
        small_molecule_ff: force field for the ligand
        protein_ff: force field for proteins
        water_ff: force field for water
        equilibration_time_ns: equilibration time per trajectory, in ns
        simulation_time_ns: simulation time per trajectory, in ns
        temperature: temperature, in K
        pressure_atm: pressure, in atm
        langevin_timescale_ps: timescale for the Langevin integrator, in ps⁻¹
        timestep_fs: timestep, in femtoseconds
        hydrogen_mass: hydrogen mass, in atomic mass units
        ligand_residue_name: name of the residue corresponding to the ligand
        constrain_hydrogens: whether to use SHAKE to freeze bonds to hydrogen
        nonbonded_cutoff: nonbonded cutoff for particle-mesh Ewald, in Å
        ionic_strength_M: ionic strength of the solution, in M (molar)
        water_buffer: amount of water to add around the protein, in Å
        protein_restraint_cutoff: cutoff past which alpha-carbons will be constrained, in Å,
            measured from the ligand. None applies no restraints
        protein_restraint_constant: force constant for backbone restraints, in kcal/mol/Å²
        save_solvent: whether to save solvent molecules
        num_solvent_to_save: number of solvent molecules to save (the N nearest the ligand each
            frame). None saves all solvent when save_solvent is True
        analysis_interval_ps: interval at which to compute per-frame SASA and polar SASA, in ps.
            None disables those analyses
        clustering: how to cluster trajectory frames. None disables clustering; pass a
            KMeansClusteringSettings (num_clusters) or GreedyClusteringSettings (cutoff_angstrom)
        validate_forcefield: validate the protein forcefield
            compatibility before submitting. Raises an error early if the protein cannot
            be parameterized or has clashing residues
        name: name of the workflow
        folder_uuid: UUID of the folder to place the workflow in
        folder: destination folder
        max_credits: maximum credits for the workflow
        webhook_url: URL that Rowan will POST to when the workflow completes
        is_draft: save as a draft without starting execution

    Returns:
        submitted workflow

    Raises:
        httpx.HTTPStatusError: request to the API fails
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
        small_molecule_ff=small_molecule_ff,
        protein_ff=protein_ff,
        water_ff=water_ff,
        num_trajectories=num_trajectories,
        equilibration_time_ns=equilibration_time_ns,
        simulation_time_ns=simulation_time_ns,
        temperature=temperature,
        pressure_atm=pressure_atm,
        langevin_timescale_ps=langevin_timescale_ps,
        timestep_fs=timestep_fs,
        hydrogen_mass=hydrogen_mass,
        ligand_residue_name=ligand_residue_name,
        constrain_hydrogens=constrain_hydrogens,
        nonbonded_cutoff=nonbonded_cutoff,
        ionic_strength_M=ionic_strength_M,
        water_buffer=water_buffer,
        protein_restraint_cutoff=protein_restraint_cutoff,
        protein_restraint_constant=protein_restraint_constant,
        save_solvent=save_solvent,
        num_solvent_to_save=num_solvent_to_save,
        analysis_interval_ps=analysis_interval_ps,
        clustering=clustering,
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
        return Workflow[PoseAnalysisMDResult](**response.json())
