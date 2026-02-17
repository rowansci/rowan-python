"""NMR workflow - predict Nuclear Magnetic Resonance spectra."""

from dataclasses import dataclass
from typing import Any

import stjames
from rdkit import Chem

from ..utils import api_client
from .base import RdkitMol, StJamesMolecule, Workflow, WorkflowResult, register_result


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
        """Per-atom NMR chemical shifts (ensemble average)."""
        return list(self._workflow.chemical_shifts)

    @property
    def boltzmann_weights(self) -> list[float]:
        """Boltzmann weights for conformers."""
        return list(self._workflow.boltzmann_weights)

    @property
    def predicted_peaks(self) -> dict[int, list[NMRPeak]]:
        """Predicted NMR peaks by nucleus atomic number."""
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
        """0-indexed atoms which are equivalent to one another."""
        return self._workflow.symmetry_equivalent_nuclei


def submit_nmr_workflow(
    initial_molecule: dict[str, Any] | StJamesMolecule | RdkitMol,
    solvent: str | None = "chloroform",
    do_csearch: bool = True,
    do_optimization: bool = True,
    name: str = "NMR Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a Nuclear Magnetic Resonance (NMR) prediction workflow to the API.

    :param initial_molecule: The molecule to predict NMR spectra for.
    :param solvent: The solvent in which to compute NMR spectra.
    :param do_csearch: Whether to perform a conformational search on the input structure.
    :param do_optimization: Whether to perform an optimization on the input structure.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to store the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if isinstance(initial_molecule, StJamesMolecule):
        initial_molecule = initial_molecule.model_dump(mode="json")
    elif isinstance(initial_molecule, Chem.rdchem.Mol | Chem.rdchem.RWMol):
        initial_molecule = StJamesMolecule.from_rdkit(initial_molecule, cid=0)

    workflow_data = {"initial_molecule": initial_molecule, "solvent": solvent}

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
        "initial_molecule": initial_molecule,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


__all__ = ["NMRPeak", "NMRResult", "submit_nmr_workflow"]
