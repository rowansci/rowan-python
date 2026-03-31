"""Protein MD workflow - molecular dynamics simulations for proteins."""

from pathlib import Path

import stjames

from ..folder import Folder
from ..protein import Protein, retrieve_protein
from ..utils import api_client
from .base import Message, Workflow, WorkflowResult, parse_messages, register_result


@register_result("protein_md")
class ProteinMDResult(WorkflowResult):
    """Result from a Protein Molecular Dynamics (MD) workflow."""

    _stjames_class = stjames.ProteinMolecularDynamicsWorkflow

    def __repr__(self) -> str:
        n_traj = len(self.trajectory_uuids)
        return f"<ProteinMDResult trajectories={n_traj}>"

    @property
    def trajectory_uuids(self) -> list[str]:
        """UUIDs of all trajectory calculations."""
        raw = getattr(self._workflow, "trajectories", []) or []
        return [t.uuid for t in raw]

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
    def bonds(self) -> list[tuple[int, int]]:
        """Bond connectivity as pairs of atom indices."""
        raw = getattr(self._workflow, "bonds", []) or []
        return [tuple(bond) for bond in raw]

    @property
    def messages(self) -> list[Message]:
        """Any messages or warnings from the workflow."""
        return parse_messages(getattr(self._workflow, "messages", None))

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


def submit_protein_md_workflow(
    protein: str | Protein,
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
    save_solvent: bool = False,
    validate_forcefield: bool = True,
    name: str = "Protein MD Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
) -> Workflow:
    """
    Submits a Protein Molecular Dynamics (MD) workflow to the API.

    :param protein: *holo* protein on which MD will be run.
        Can be input as a UUID or a Protein object.
    :param num_trajectories: Number of trajectories to run.
    :param equilibration_time_ns: how long to equilibrate trajectories for, in ns
    :param simulation_time_ns: how long to run trajectories for, in ns
    :param temperature: temperature, in K
    :param pressure_atm: pressure, in atm
    :param langevin_timescale_ps: timescale for the Langevin integrator, in ps^-1
    :param timestep_fs: timestep, in femtoseconds
    :param constrain_hydrogens: whether or not to use SHAKE to freeze bonds to hydrogen
    :param nonbonded_cutoff: nonbonded cutoff for particle-mesh Ewald, in A
    :param ionic_strength_M: ionic strength of the solution, in M (molar)
    :param water_buffer: amount of water to add around the protein, in A
    :param save_solvent: whether solvent should be saved
    :param validate_forcefield: if True (default), validate the protein forcefield
        compatibility before submitting. Raises an error early if the protein cannot
        be parameterized or has clashing residues.
    :param name: Name of the workflow.
    :param folder_uuid: UUID of the folder to place the workflow in.
    :param folder: Folder object to store the workflow in.
    :param max_credits: Maximum number of credits to use for the workflow.
    :param webhook_url: URL that Rowan will POST to when the workflow completes.
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

    workflow = stjames.ProteinMolecularDynamicsWorkflow(
        protein=protein,
        num_trajectories=num_trajectories,
        equilibration_time_ns=equilibration_time_ns,
        simulation_time_ns=simulation_time_ns,
        temperature=temperature,
        pressure_atm=pressure_atm,
        langevin_timescale_ps=langevin_timescale_ps,
        timestep_fs=timestep_fs,
        constrain_hydrogens=constrain_hydrogens,
        nonbonded_cutoff=nonbonded_cutoff,
        ionic_strength_M=ionic_strength_M,
        water_buffer=water_buffer,
        save_solvent=save_solvent,
    )

    data = {
        "workflow_type": "protein_md",
        "workflow_data": workflow.model_dump(serialize_as_any=True, mode="json"),
        "name": name,
        "folder_uuid": folder_uuid,
        "max_credits": max_credits,
        "webhook_url": webhook_url,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())
