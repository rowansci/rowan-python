"""Pose-analysis MD workflow - molecular dynamics simulations for ligand-protein complexes."""

import stjames

from ..protein import Protein
from ..utils import api_client
from .base import Workflow, WorkflowResult, register_result


@register_result("pose_analysis_md")
class PoseAnalysisMDResult(WorkflowResult):
    """Result from a Pose-Analysis Molecular Dynamics (MD) workflow."""

    _stjames_class = stjames.PoseAnalysisMolecularDynamicsWorkflow

    def __repr__(self) -> str:
        n_traj = getattr(self._workflow, "num_trajectories", None)
        sim_time = getattr(self._workflow, "simulation_time_ns", None)
        return f"<PoseAnalysisMDResult trajectories={n_traj} simulation_ns={sim_time}>"


def submit_pose_analysis_md_workflow(
    protein: str | Protein,
    initial_smiles: str,
    num_trajectories: int = 1,
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
    protein_restraint_cutoff: float | None = None,
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


__all__ = ["PoseAnalysisMDResult", "submit_pose_analysis_md_workflow"]
