"""Shared result behavior for molecular-dynamics workflows."""

import struct
from collections.abc import Sequence
from pathlib import Path
from typing import Protocol

from stjames.pdb import PDB, pdb_object_to_pdb_filestring

from ..protein import Protein, retrieve_protein
from ..utils import api_client, download_file
from .base import WorkflowResult


class _MolecularDynamicsTrajectory(Protocol):
    """Trajectory fields used by shared molecular-dynamics result methods."""

    @property
    def mean_structure_uuid(self) -> str | None:
        """UUID of the coordinate-averaged structure."""
        ...

    @property
    def median_structure_frame_index(self) -> int | None:
        """Frame index of the medoid structure."""
        ...


class _MolecularDynamicsResult(WorkflowResult):
    """Base result with shared molecular-dynamics retrieval and download methods."""

    @property
    def trajectories(self) -> Sequence[_MolecularDynamicsTrajectory]:
        """Per-replicate trajectory results."""
        raise NotImplementedError

    @property
    def minimized_protein_uuid(self) -> str | None:
        """UUID of the energy-minimized protein structure."""
        return getattr(self._workflow, "minimized_protein_uuid", None)

    def get_minimized_protein(self) -> Protein | None:
        """Fetch and cache the energy-minimized protein structure.

        :returns: protein object, or `None` when no minimized structure is available
        """
        if not (uuid := self.minimized_protein_uuid):
            return None
        if "minimized_protein" not in self._cache:
            self._cache["minimized_protein"] = retrieve_protein(
                uuid, workflow_uuid=self.workflow_uuid
            )
        return self._cache["minimized_protein"]

    def get_mean_structure(self, replicate: int = 0) -> Protein | None:
        """Fetch and cache the coordinate-averaged structure for a replicate.

        :param replicate: zero-based trajectory replicate index
        :returns: protein object, or `None` when no mean structure is available
        :raises IndexError: if `replicate` is out of range
        """
        trajectory = self._trajectory(replicate)
        if not trajectory.mean_structure_uuid:
            return None
        cache_key = f"mean_structure_{replicate}"
        if cache_key not in self._cache:
            self._cache[cache_key] = retrieve_protein(
                trajectory.mean_structure_uuid,
                workflow_uuid=self.workflow_uuid,
            )
        return self._cache[cache_key]

    def download_mean_structure(
        self,
        replicate: int = 0,
        path: Path | str | None = None,
        name: str | None = None,
    ) -> Path | None:
        """Download the coordinate-averaged structure for a replicate as PDB.

        :param replicate: zero-based trajectory replicate index
        :param path: output directory; defaults to the current directory
        :param name: filename without the `.pdb` extension
        :returns: downloaded path, or `None` when no mean structure is available
        :raises IndexError: if `replicate` is out of range
        """
        protein = self.get_mean_structure(replicate)
        if protein is None:
            return None
        directory = Path(path) if path is not None else Path.cwd()
        file_name = name or f"mean_structure_{replicate}"
        protein.download_pdb_file(
            path=directory,
            name=file_name,
            workflow_uuid=self.workflow_uuid,
        )
        return directory / f"{file_name}.pdb"

    def download_medoid_structure(
        self,
        replicate: int = 0,
        path: Path | str | None = None,
        name: str | None = None,
    ) -> Path | None:
        """Download the medoid trajectory frame for a replicate as PDB.

        :param replicate: zero-based trajectory replicate index
        :param path: output directory; defaults to the current directory
        :param name: filename without the `.pdb` extension
        :returns: downloaded path, or `None` when no medoid frame is available
        :raises IndexError: if `replicate` is out of range
        :raises ValueError: if topology and trajectory data are inconsistent
        """
        trajectory = self._trajectory(replicate)
        if trajectory.median_structure_frame_index is None:
            return None
        minimized = self.get_minimized_protein()
        if minimized is None:
            return None
        if minimized.data is None:
            minimized.refresh(workflow_uuid=self.workflow_uuid)
        if minimized.data is None:
            raise ValueError("Minimized protein response did not include structure data.")

        pdb_string = pdb_object_to_pdb_filestring(
            pdb=PDB.model_validate(minimized.data),
            header=True,
            source=True,
            keyword=True,
            crystallography=True,
            remark=False,
        )
        atom_count = sum(line.startswith(("ATOM  ", "HETATM")) for line in pdb_string.splitlines())
        with api_client() as client:
            response = client.get(
                f"/trajectory/{self.workflow_uuid}/compressed_stream",
                params={
                    "replicate": replicate,
                    "start_frame": trajectory.median_structure_frame_index,
                    "num_frames": 1,
                },
            )
            response.raise_for_status()

        box_size = int(response.headers.get("X-Box-Size-Bytes", 48))
        response_atoms = int(response.headers.get("X-Num-Atoms", atom_count))
        if response_atoms != atom_count:
            raise ValueError(
                f"Trajectory has {response_atoms} atoms but minimized topology has {atom_count}."
            )
        coordinate_count = atom_count * 3
        coordinate_bytes = response.content[box_size : box_size + coordinate_count * 8]
        if len(coordinate_bytes) != coordinate_count * 8:
            raise ValueError("Trajectory response did not contain one complete coordinate frame.")
        coordinates = iter(struct.unpack(f"<{coordinate_count}d", coordinate_bytes))

        output_lines: list[str] = []
        for line in pdb_string.splitlines():
            if line.startswith(("ATOM  ", "HETATM")):
                x, y, z = next(coordinates), next(coordinates), next(coordinates)
                line = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
            output_lines.append(line)

        directory = Path(path) if path is not None else Path.cwd()
        directory.mkdir(parents=True, exist_ok=True)
        file_path = directory / f"{name or f'medoid_structure_{replicate}'}.pdb"
        file_path.write_text("\n".join(output_lines) + "\n")
        return file_path

    def get_atom_distances(
        self,
        atom_pairs: list[tuple[int, int]],
        replicate: int = 0,
    ) -> list[list[float]]:
        """Fetch interatomic distances over a trajectory.

        :param atom_pairs: zero-based atom-index pairs
        :param replicate: zero-based trajectory replicate index
        :returns: distance arrays in angstrom, one per atom pair
        :raises requests.HTTPError: if the API request fails
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
        """Download DCD trajectory files for selected replicates.

        :param replicates: zero-based replicate indices
        :param name: archive filename without the `.tar.gz` extension
        :param path: output directory; defaults to the current directory
        :returns: path to the downloaded `.tar.gz` archive
        :raises requests.HTTPError: if the API request fails
        """
        directory = Path(path) if path is not None else Path.cwd()
        file_path = directory / f"{name or 'trajectories'}.tar.gz"
        return download_file(
            file_path,
            "POST",
            f"/trajectory/{self.workflow_uuid}/trajectory_dcds",
            json=replicates,
        )

    def _trajectory(self, replicate: int) -> _MolecularDynamicsTrajectory:
        """Return one trajectory with a consistent bounds error."""
        trajectories = self.trajectories
        if replicate < 0 or replicate >= len(trajectories):
            raise IndexError(f"Replicate {replicate} out of range (0-{len(trajectories) - 1})")
        return trajectories[replicate]
