"""Tautomer-search workflow - find tautomeric forms of molecules."""

from dataclasses import dataclass
from typing import Any

import stjames

from ..calculation import retrieve_calculation
from ..folder import Folder
from ..molecule import Molecule
from ..utils import api_client
from .base import (
    Message,
    Mode,
    MoleculeInput,
    Workflow,
    WorkflowResult,
    molecule_to_dict,
    parse_messages,
    register_result,
)


@dataclass(frozen=True, slots=True)
class Tautomer:
    """
    Tautomer result.

    :param energy: Energy in Hartree.
    :param weight: Boltzmann weight (sum to 1.0 across all tautomers).
    :param predicted_relative_energy: Relative energy in kcal/mol (relative to lowest energy).
    :param structure_uuids: UUIDs of the structure calculations.
    """

    energy: float
    weight: float
    predicted_relative_energy: float
    structure_uuids: tuple[str, ...]


@register_result("tautomers")
class TautomerResult(WorkflowResult):
    """Result from a tautomer-search workflow."""

    _stjames_class = stjames.TautomerWorkflow

    def __post_init__(self) -> None:
        super().__post_init__()
        tautomers = self.tautomers
        if tautomers:
            best = max(tautomers, key=lambda t: t.weight)
            if best.structure_uuids:
                calc = retrieve_calculation(best.structure_uuids[0])
                self._cache["best_tautomer"] = calc.molecule

    def __repr__(self) -> str:
        tautomers = self.tautomers
        n = len(tautomers)
        if tautomers:
            lowest = min(tautomers, key=lambda t: t.energy)
            return f"<TautomerResult tautomers={n} lowest_energy={lowest.energy:.6f} H>"
        return f"<TautomerResult tautomers={n}>"

    @property
    def tautomers(self) -> list[Tautomer]:
        """List of tautomers with energies and weights."""
        return [
            Tautomer(
                energy=t.energy,
                weight=t.weight,
                predicted_relative_energy=t.predicted_relative_energy,
                structure_uuids=tuple(s.uuid for s in t.structures if s.uuid),
            )
            for t in self._workflow.tautomers
        ]

    @property
    def messages(self) -> list[Message]:
        """Any messages or warnings from the workflow."""
        return parse_messages(getattr(self._workflow, "messages", None))

    @property
    def best_tautomer(self) -> Molecule | None:
        """Molecule of the highest Boltzmann-weight tautomer."""
        return self._cache.get("best_tautomer")

    @property
    def molecules(self) -> list[Molecule]:
        """Molecules for all tautomers.

        .. note::
            Makes one API call per tautomer on first access.
            Results are cached. Call clear_cache() to refresh.
        """
        if "all_molecules" not in self._cache:
            mols = []
            for t in self.tautomers:
                if t.structure_uuids:
                    calc = retrieve_calculation(t.structure_uuids[0])
                    if calc.molecule:
                        mols.append(calc.molecule)
            self._cache["all_molecules"] = mols
        return self._cache["all_molecules"]


def submit_tautomer_search_workflow(
    initial_molecule: MoleculeInput,
    mode: Mode = Mode.CAREFUL,
    conf_gen_settings: stjames.ConformerGenSettingsUnion | None = None,
    multistage_opt_settings: stjames.MultiStageOptSettings | None = None,
    name: str = "Tautomer Search Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
    is_draft: bool = False,
) -> Workflow:
    """
    Submits a tautomer-search workflow to the API.

    :param initial_molecule: Molecule to find tautomers for.
    :param mode: *Deprecated.* Ignored when ``multistage_opt_settings`` is provided.
    :param conf_gen_settings: Conformer generation settings. Defaults to ETKDG with
        250 initial conformers, 20 max conformers, and 15 kcal/mol MMFF energy cutoff.
    :param multistage_opt_settings: Optimization settings for tautomer ranking.
        Defaults to AIMNet2/wB97M-D3 optimization with CPCMx singlepoint.
    :param name: Name of the workflow.
    :param folder_uuid: UUID of the folder to place the workflow in.
    :param folder: Folder object to store the workflow in.
    :param max_credits: Maximum number of credits to use for the workflow.
    :param webhook_url: URL that Rowan will POST to when the workflow completes.
    :param is_draft: If True, submit the workflow as a draft without starting execution.
    :returns: Workflow object representing the submitted workflow.
    :raises requests.HTTPError: If the request to the API fails.
    """
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder:
        folder_uuid = folder.uuid
    initial_molecule = molecule_to_dict(initial_molecule)

    workflow_kwargs: dict[str, Any] = {"initial_molecule": initial_molecule, "mode": mode}
    if conf_gen_settings is not None:
        workflow_kwargs["conf_gen_settings"] = conf_gen_settings
    if multistage_opt_settings is not None:
        workflow_kwargs["multistage_opt_settings"] = multistage_opt_settings

    workflow = stjames.TautomerWorkflow(**workflow_kwargs)

    data = {
        "workflow_type": "tautomers",
        "workflow_data": workflow.model_dump(serialize_as_any=True, mode="json"),
        "initial_molecule": initial_molecule,
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
