"""NMR workflow - predict Nuclear Magnetic Resonance spectra."""

from dataclasses import dataclass
from typing import Any

import stjames

from ..folder import Folder
from ..types import SolventInput
from ..utils import api_client
from .base import (
    StructureInput,
    Workflow,
    WorkflowResult,
    molecule_to_dict,
    register_result,
    require_coordinates,
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


def _nmr_multistage_opt_settings(solvent: SolventInput) -> stjames.MultiStageOptSettings:
    """
    Build NMR optimization settings, adding a solvated AIMNet2 singlepoint.

    Optimization runs gas-phase (matching MagNet); the solvated singlepoint, whose solvent model
    is looked up from stjames.NMR_SOLVENT_MODELS, reweights the conformer ensemble.

    :param solvent: solvent for the prediction
    :returns: multi-stage optimization settings
    """
    singlepoint_settings = stjames.Settings(
        method=stjames.Method.AIMNET2_WB97MD3,
        tasks=[stjames.Task.ENERGY],
        solvent_settings=stjames.SolventSettings(
            solvent=solvent, model=stjames.NMR_SOLVENT_MODELS[stjames.Solvent(solvent)]
        ),
    )

    return stjames.MultiStageOptSettings(
        optimization_settings=[
            stjames.Settings(method=stjames.Method.AIMNET2_WB97MD3, tasks=[stjames.Task.OPTIMIZE])
        ],
        singlepoint_settings=singlepoint_settings,
    )


def submit_nmr_workflow(
    initial_molecule: StructureInput,
    solvent: SolventInput = "chloroform",
    do_csearch: bool = False,
    do_optimization: bool = True,
    name: str = "NMR Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
    is_draft: bool = False,
) -> Workflow:
    """
    Submits a Nuclear Magnetic Resonance (NMR) prediction workflow to the API.

    :param initial_molecule: Molecule to predict NMR spectra for.
    :param solvent: Solvent for NMR calculation (default: chloroform). Must be an NMR-supported
        solvent (see rowan.NMR_SUPPORTED_SOLVENTS); others raise ValueError. A solvated AIMNet2
        singlepoint reweights the conformer ensemble, using CPCM-X where supported and
        otherwise ALPB.
    :param do_csearch: Whether to perform a conformational search. Requires do_optimization.
    :param do_optimization: Whether to optimize conformer geometries.
    :param name: Name of the workflow.
    :param folder_uuid: UUID of the folder to store the workflow in.
    :param folder: Folder object to store the workflow in.
    :param max_credits: Maximum number of credits to use for the workflow.
    :param webhook_url: URL that Rowan will POST to when the workflow completes.
    :param is_draft: If True, submit the workflow as a draft without starting execution.
    :returns: Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    require_coordinates(initial_molecule)
    if stjames.Solvent(solvent) not in stjames.NMR_SUPPORTED_SOLVENTS:
        supported = ", ".join(sorted(s.value for s in stjames.NMR_SUPPORTED_SOLVENTS))
        raise ValueError(
            f"{stjames.Solvent(solvent).value!r} is not an NMR-supported solvent. "
            f"NMR-supported solvents: {supported}."
        )
    if do_csearch and not do_optimization:
        raise ValueError(
            "`do_optimization` must be True when `do_csearch` is True; the conformers from "
            "the search must be optimized before NMR prediction."
        )
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder:
        folder_uuid = folder.uuid
    mol_dict = molecule_to_dict(initial_molecule)

    workflow_data: dict[str, Any] = {"initial_molecule": mol_dict, "solvent": solvent}

    if not do_csearch:
        workflow_data["conf_gen_settings"] = None

    if not do_optimization:
        workflow_data["multistage_opt_settings"] = None
    else:
        workflow_data["multistage_opt_settings"] = _nmr_multistage_opt_settings(solvent)

    workflow = stjames.NMRSpectroscopyWorkflow.model_validate(workflow_data)

    data = {
        "workflow_type": "nmr",
        "workflow_data": workflow.model_dump(serialize_as_any=True, mode="json"),
        "initial_molecule": mol_dict,
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
