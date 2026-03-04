"""MSA workflow - Multiple Sequence Alignment for proteins."""

from pathlib import Path
from typing import Literal

import stjames

from ..utils import api_client
from .base import Workflow, WorkflowResult, register_result

MSAOutputFormat = Literal["colabfold", "chai", "boltz"]


@register_result("msa")
class MSAResult(WorkflowResult):
    """Result from a Multiple Sequence Alignment (MSA) workflow."""

    _stjames_class = stjames.MSAWorkflow

    def __repr__(self) -> str:
        return f"<MSAResult uuid='{self.workflow_uuid}'>"

    @property
    def _output_formats(self) -> list[str]:
        """Output formats requested for the MSA (internal use)."""
        fmts = getattr(self._workflow, "output_formats", []) or []
        return [f.value if hasattr(f, "value") else str(f) for f in fmts]

    def download_files(
        self,
        format: MSAOutputFormat | None = None,
        path: Path | None = None,
    ) -> list[Path]:
        """
        Download MSA files for this workflow.

        :param format: Output format to download. If None, downloads all requested formats.
        :param path: Directory to save files to. Defaults to current directory.
        :return: List of paths to downloaded tar.gz files.
        :raises ValueError: If the requested format wasn't in the original output_formats.
        :raises HTTPError: If the API request fails.
        """
        if format is not None and format not in self._output_formats:
            raise ValueError(
                f"Format '{format}' was not requested. Available formats: {self._output_formats}"
            )

        if path is None:
            path = Path.cwd()

        path.mkdir(parents=True, exist_ok=True)

        formats_to_download = [format] if format else self._output_formats
        downloaded_paths = []

        for fmt in formats_to_download:
            with api_client() as client:
                response = client.get(
                    f"/workflow/{self.workflow_uuid}/get_msa_files",
                    params={"msa_format": fmt},
                )
                response.raise_for_status()

            file_path = path / f"msa-{fmt}.tar.gz"
            with open(file_path, "wb") as f:
                f.write(response.content)
            downloaded_paths.append(file_path)

        return downloaded_paths


def submit_msa_workflow(
    initial_protein_sequences: list[str | stjames.ProteinSequence],
    output_formats: list[MSAOutputFormat | stjames.MSAFormat] | None = None,
    name: str = "MSA Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a Multiple Sequence Alignment (MSA) workflow to the API.

    :param initial_protein_sequences: List of protein sequences to align (amino acid strings).
    :param output_formats: Output formats for the MSA files. Options: "colabfold", "chai", "boltz".
    :param name: The name to assign to the workflow.
    :param folder_uuid: UUID of the folder where the workflow will be stored.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted MSA workflow.
    :raises HTTPError: If the API request fails.
    """
    if output_formats is None:
        output_formats = ["colabfold"]

    # Convert to stjames types
    protein_sequences = []
    for seq in initial_protein_sequences:
        if isinstance(seq, stjames.ProteinSequence):
            protein_sequences.append(seq)
        else:
            protein_sequences.append(stjames.ProteinSequence(sequence=seq))

    msa_formats = []
    for fmt in output_formats:
        if isinstance(fmt, stjames.MSAFormat):
            msa_formats.append(fmt)
        else:
            msa_formats.append(stjames.MSAFormat(fmt))

    workflow = stjames.MSAWorkflow(
        initial_protein_sequences=protein_sequences,
        output_formats=msa_formats,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "msa",
        "workflow_data": workflow.model_dump(serialize_as_any=True, mode="json"),
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


__all__ = ["MSAOutputFormat", "MSAResult", "submit_msa_workflow"]
