"""NMR workflow - predict Nuclear Magnetic Resonance spectra."""

from dataclasses import dataclass

import stjames

from ..types import MoleculeInput, SolventInput
from ..utils import api_client
from .base import (
    Workflow,
    WorkflowResult,
    molecule_to_dict,
    register_result,
)


@dataclass(frozen=True, slots=True)
class NMRPeak:
    """An NMR peak."""

    nucleus: int
    shift: float
    atom_indices: tuple[int, ...]


@register_result("nmr")
class NMRResult(WorkflowResult):
    """Result from a Nuclear Magnetic Resonance (NMR) workflow."""

    _stjames_class = stjames.NMRSpectroscopyWorkflow

    def __repr__(self) -> str:
        shifts = self.chemical_shifts
        n_atoms = len(shifts)
        n_conformers = len(self.boltzmann_weights)
        return f"<NMRResult atoms={n_atoms} conformers={n_conformers}>"

    @property
    def chemical_shifts(self) -> list[float | None]:
        """
        Per-atom NMR chemical shifts (Boltzmann-weighted ensemble average).

        Index corresponds to atom index in the molecule. Returns None for
        atoms without NMR-active nuclei (e.g., oxygen).
        """
        return list(self._workflow.chemical_shifts)

    @property
    def per_conformer_chemical_shifts(self) -> list[list[float | None]]:
        """
        Chemical shifts for each conformer before Boltzmann averaging.

        Outer list is per-conformer, inner list is per-atom.
        """
        return [list(shifts) for shifts in self._workflow.per_conformer_chemical_shifts]

    @property
    def boltzmann_weights(self) -> list[float]:
        """Boltzmann weights for each conformer (sum to 1.0)."""
        return list(self._workflow.boltzmann_weights)

    @property
    def conformer_uuids(self) -> list[str]:
        """UUIDs of the conformer calculations."""
        return self._workflow.conformers

    @property
    def predicted_peaks(self) -> dict[int, list[NMRPeak]]:
        """
        Predicted NMR peaks grouped by nucleus atomic number.

        Keys are atomic numbers (1 for 1H, 6 for 13C). Peaks with equivalent
        atoms are merged and shifts are averaged.
        """
        return {
            nucleus: [
                NMRPeak(
                    nucleus=p.nucleus,
                    shift=p.shift,
                    atom_indices=tuple(p.atom_indices),
                )
                for p in peaks
            ]
            for nucleus, peaks in self._workflow.predicted_peaks.items()
        }

    @property
    def symmetry_equivalent_nuclei(self) -> list[list[int]]:
        """
        Groups of symmetry-equivalent atom indices (0-indexed).

        Atoms in the same group have equivalent chemical environments
        and are averaged together in predicted_peaks.
        """
        return self._workflow.symmetry_equivalent_nuclei


def submit_nmr_workflow(
    initial_molecule: MoleculeInput,
    solvent: SolventInput = "chloroform",
    do_csearch: bool = True,
    do_optimization: bool = True,
    name: str = "NMR Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a Nuclear Magnetic Resonance (NMR) prediction workflow to the API.

    :param initial_molecule: Molecule to predict NMR spectra for.
    :param solvent: Solvent for NMR calculation (default: chloroform).
    :param do_csearch: Whether to perform a conformational search.
    :param do_optimization: Whether to optimize conformer geometries.
    :param name: Name of the workflow.
    :param folder_uuid: UUID of the folder to store the workflow in.
    :param max_credits: Maximum number of credits to use for the workflow.
    :returns: Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    mol_dict = molecule_to_dict(initial_molecule)

    workflow_data = {"initial_molecule": mol_dict, "solvent": solvent}

    if not do_csearch:
        workflow_data["conf_gen_settings"] = None

    if not do_optimization:
        workflow_data["multistage_opt_settings"] = None

    workflow = stjames.NMRSpectroscopyWorkflow.model_validate(workflow_data)

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "nmr",
        "workflow_data": workflow.model_dump(serialize_as_any=True, mode="json"),
        "initial_molecule": mol_dict,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())
