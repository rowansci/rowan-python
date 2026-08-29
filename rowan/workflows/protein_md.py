"""Protein MD workflow - molecular dynamics simulations for proteins."""

from dataclasses import dataclass, field
from typing import Literal

import stjames
from stjames import (
    Binder,
    GreedyClusteringSettings,
    KMeansClusteringSettings,
    ProteinForceField,
    WaterForceField,
)

from ..folder import Folder
from ..protein import Protein
from ..types import ProteinUUID
from ..utils import api_client
from ._molecular_dynamics import _MolecularDynamicsResult
from .base import Message, Workflow, parse_messages, register_result


@dataclass(frozen=True, slots=True)
class ProteinMDTrajectory:
    """
    Results from a single protein-MD trajectory replicate.

    :param uuid: UUID of the trajectory calculation.
    :param sasa: Solvent-accessible surface area per analyzed frame (populated when
        analysis_interval_ps is set).
    :param polar_sasa: Polar solvent-accessible surface area per analyzed frame (populated when
        analysis_interval_ps is set).
    :param isotropic_radius_of_gyration: Radius of gyration per analyzed frame.
    :param cluster_centroid_indices: Frame indices of the cluster centroids (populated when
        clustering is set).
    :param cluster_indices_by_frame: Cluster assignment for each frame (populated when clustering
        is set).
    :param binder_rmsd: Per-frame binder RMSD, when the binder has one component.
    :param mmgbsa_scores: Per-frame MM/GBSA score for the complete binder.
    :param mean_structure_uuid: UUID of the coordinate-averaged structure.
    :param median_structure_frame_index: Frame index of the medoid structure.
    """

    uuid: str
    sasa: list[float | None]
    polar_sasa: list[float | None]
    isotropic_radius_of_gyration: list[float]
    cluster_centroid_indices: list[int]
    cluster_indices_by_frame: list[int]
    binder_rmsd: list[float] = field(default_factory=list)
    mmgbsa_scores: list[float | None] = field(default_factory=list)
    mean_structure_uuid: str | None = None
    median_structure_frame_index: int | None = None


@register_result("protein_md")
class ProteinMDResult(_MolecularDynamicsResult):
    """Result from a Protein Molecular Dynamics (MD) workflow."""

    _stjames_class = stjames.ProteinMolecularDynamicsWorkflow

    def __post_init__(self) -> None:
        """Normalize the pre-0.0.255 binder schema before parsing old workflows."""
        binder = self.workflow_data.get("binder")
        if isinstance(binder, dict) and "small_molecules" in binder:
            small_molecules = binder.pop("small_molecules")
            if isinstance(small_molecules, dict):
                self.workflow_data.setdefault("small_molecules", small_molecules)
                binder["small_molecule_residues"] = list(small_molecules)
            elif isinstance(small_molecules, list):
                binder["small_molecule_residues"] = small_molecules
        super().__post_init__()

    def __repr__(self) -> str:
        n_traj = len(self.trajectory_uuids)
        return f"<ProteinMDResult trajectories={n_traj}>"

    @property
    def trajectory_uuids(self) -> list[str]:
        """UUIDs of all trajectory calculations."""
        raw = getattr(self._workflow, "trajectories", []) or []
        return [t.uuid for t in raw]

    @property
    def trajectories(self) -> list[ProteinMDTrajectory]:
        """Per-replicate trajectory results (SASA, radius of gyration, cluster assignments)."""
        raw = getattr(self._workflow, "trajectories", []) or []
        return [
            ProteinMDTrajectory(
                uuid=t.uuid,
                sasa=t.sasa,
                polar_sasa=t.polar_sasa,
                isotropic_radius_of_gyration=t.isotropic_radius_of_gyration,
                cluster_centroid_indices=t.cluster_centroid_indices,
                cluster_indices_by_frame=t.cluster_indices_by_frame,
                binder_rmsd=t.binder_rmsd,
                mmgbsa_scores=t.mmgbsa_scores,
                mean_structure_uuid=t.mean_structure_uuid,
                median_structure_frame_index=t.median_structure_frame_index,
            )
            for t in raw
        ]

    @property
    def bonds(self) -> list[tuple[int, int]]:
        """Bond connectivity as pairs of atom indices."""
        raw = getattr(self._workflow, "bonds", []) or []
        return [tuple(bond) for bond in raw]

    @property
    def messages(self) -> list[Message]:
        """Any messages or warnings from the workflow."""
        return parse_messages(getattr(self._workflow, "messages", None))


def submit_protein_md_workflow(
    protein: Protein | ProteinUUID,
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
    save_solvent: bool = False,
    num_solvent_to_save: int | None = None,
    small_molecules: dict[str | int, str | None] | None = None,
    binder: Binder | None = None,
    protein_restraint_cutoff: float | None = None,
    protein_restraint_constant: float = 100,
    analysis_interval_ps: float | None = None,
    clustering: KMeansClusteringSettings | GreedyClusteringSettings | None = None,
    validate_forcefield: bool = True,
    name: str = "Protein MD Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
    is_draft: bool = False,
) -> Workflow:
    """
    Submits a Protein Molecular Dynamics (MD) workflow to the API.

    :param protein: *holo* protein on which MD will be run.
        Can be input as a UUID or a Protein object.
    :param num_trajectories: Number of trajectories to run.
    :param small_molecule_ff: Force field for small molecules.
    :param protein_ff: Force field for proteins.
    :param water_ff: Force field for water.
    :param equilibration_time_ns: how long to equilibrate trajectories for, in ns
    :param simulation_time_ns: how long to run trajectories for, in ns
    :param temperature: temperature, in K
    :param pressure_atm: pressure, in atm
    :param langevin_timescale_ps: timescale for the Langevin integrator, in ps^-1
    :param timestep_fs: timestep, in femtoseconds
    :param hydrogen_mass: hydrogen mass, in atomic mass units
    :param constrain_hydrogens: whether or not to use SHAKE to freeze bonds to hydrogen
    :param nonbonded_cutoff: nonbonded cutoff for particle-mesh Ewald, in A
    :param ionic_strength_M: ionic strength of the solution, in M (molar)
    :param water_buffer: amount of water to add around the protein, in A
    :param save_solvent: whether solvent should be saved
    :param num_solvent_to_save: keep this many solvent molecules nearest the binder, or all if None;
        only meaningful when save_solvent is True and a binder is present
    :param small_molecules: SMILES by protein residue name or index for small molecules that
        require separate parameterization. A None value uses an existing residue template.
    :param binder: optional binder specification (protein chains and/or small molecules).
        When set, per-frame MM/GBSA scores are computed against the whole binder.
        Per-frame binder RMSD is populated only when the binder is a single component
        (one small molecule → heavy-atom RMSD; one binder chain → backbone N/CA/C/O RMSD);
        it is empty for multi-molecule, multi-chain, or combined chain+molecule binders.
    :param protein_restraint_cutoff: cutoff distance from the binder past which Cα atoms are
        harmonically restrained, in Å; None disables restraints
    :param protein_restraint_constant: force constant for Cα backbone restraints, in kcal/mol/Å²
    :param analysis_interval_ps: Interval at which to compute per-frame SASA and polar SASA, in ps.
        None disables those analyses.
    :param clustering: How to cluster trajectory frames. None disables clustering; pass a
        KMeansClusteringSettings (num_clusters) or GreedyClusteringSettings (cutoff_angstrom).
    :param validate_forcefield: if True (default), validate the protein forcefield
        compatibility before submitting. Raises an error early if the protein cannot
        be parameterized or has clashing residues. Binder small molecules are skipped,
        whether given by residue name or by index, since they are parameterized from their
        SMILES rather than the protein forcefield; cofactors, metals, and glycans outside
        the binder are still validated.
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
        exclude_residues = (
            list(binder.small_molecule_residues)
            if binder is not None and binder.small_molecule_residues
            else None
        )
        Protein(uuid=protein).validate_protein_forcefield(exclude_residues=exclude_residues)

    workflow = stjames.ProteinMolecularDynamicsWorkflow(
        protein=protein,
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
        constrain_hydrogens=constrain_hydrogens,
        nonbonded_cutoff=nonbonded_cutoff,
        ionic_strength_M=ionic_strength_M,
        water_buffer=water_buffer,
        save_solvent=save_solvent,
        num_solvent_to_save=num_solvent_to_save,
        small_molecules=small_molecules,
        binder=binder,
        protein_restraint_cutoff=protein_restraint_cutoff,
        protein_restraint_constant=protein_restraint_constant,
        analysis_interval_ps=analysis_interval_ps,
        clustering=clustering,
    )

    data = {
        "workflow_type": "protein_md",
        "workflow_data": workflow.model_dump(serialize_as_any=True, mode="json"),
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
