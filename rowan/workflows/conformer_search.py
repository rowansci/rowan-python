"""Conformer-search workflow - find low-energy molecular conformations."""

import stjames

from ..calculation import Calculation, retrieve_calculation
from ..folder import Folder
from ..molecule import Molecule
from ..types import MoleculeInput, SolventInput
from ..utils import api_client
from .base import (
    Workflow,
    WorkflowResult,
    molecule_to_dict,
    register_result,
)
from .constants import to_relative_kcal


@register_result("conformer_search")
class ConformerSearchResult(WorkflowResult):
    """Result from a conformer-search workflow."""

    _stjames_class = stjames.ConformerSearchWorkflow

    def __repr__(self) -> str:
        n = self.num_conformers
        return f"<ConformerSearchResult conformers={n}>"

    @property
    def num_conformers(self) -> int:
        """Number of conformers found."""
        return len(self._workflow.conformer_uuids)

    @property
    def conformer_uuids(self) -> list[list[str | None]]:
        """List of conformer UUIDs (nested for multistage optimization)."""
        return self._workflow.conformer_uuids

    def get_energies(self, relative: bool = False) -> list[float]:
        """
        Get conformer energies.

        :param relative: If True, return relative energies in kcal/mol (relative to
            the lowest energy conformer). If False (default), return absolute
            energies in Hartree.
        :returns: List of conformer energies ordered by energy (lowest first).
        """
        energies: list[float] = list(self._workflow.energies)
        return to_relative_kcal(energies) if relative else energies

    @property
    def radii_of_gyration(self) -> list[float]:
        """Radius of gyration for each conformer (A)."""
        return [p.radius_of_gyration for p in self._workflow.conformer_properties]

    @property
    def sasa(self) -> list[float]:
        """Solvent accessible surface area for each conformer (A^2)."""
        return [p.solvent_accessible_surface_area for p in self._workflow.conformer_properties]

    @property
    def polar_sasa(self) -> list[float]:
        """Polar solvent accessible surface area for each conformer (A^2)."""
        return [
            p.polar_solvent_accessible_surface_area for p in self._workflow.conformer_properties
        ]

    def get_conformers(self, n: int | None = None) -> list[Molecule]:
        """
        Fetch conformer molecules.

        :param n: Number of conformers to fetch (default: all). Conformers are
            ordered by energy, so n=5 returns the 5 lowest-energy conformers.
        :returns: List of Molecule objects.

        .. note::
            Makes one API call per conformer.
        """
        uuids = self._workflow.conformer_uuids
        count = len(uuids) if n is None else min(n, len(uuids))
        molecules = []
        for i in range(count):
            calc = self.get_conformer(i)
            if calc.molecule:
                molecules.append(calc.molecule)
        return molecules

    def get_conformer(self, index: int, stage: int = -1) -> Calculation:
        """
        Fetch a conformer's calculation data by index.

        .. note::
            Makes one API call per conformer on first access.
            Results are cached. Call clear_cache() to refresh.

        :param index: Conformer index (0-based).
        :param stage: Optimization stage (-1 for final stage).
        :returns: Calculation object with molecule and energy data.
        :raises IndexError: If the index is out of range.
        :raises ValueError: If the conformer UUID is None.
        """
        uuids = self._workflow.conformer_uuids
        if index < 0 or index >= len(uuids):
            raise IndexError(f"Conformer index {index} out of range (0-{len(uuids) - 1})")

        stage_uuids = uuids[index]
        uuid = stage_uuids[stage]
        if uuid is None:
            raise ValueError(f"Conformer {index} has no calculation at stage {stage}")

        cache_key = f"conformer_{index}_{stage}"
        if cache_key not in self._cache:
            self._cache[cache_key] = retrieve_calculation(uuid)
        return self._cache[cache_key]


def submit_conformer_search_workflow(
    initial_molecule: MoleculeInput,
    conf_gen_settings: stjames.ConformerGenSettings | None = None,
    final_method: stjames.Method | str = "aimnet2_wb97md3",
    solvent: SolventInput = None,
    transition_state: bool = False,
    name: str = "Conformer Search Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
    is_draft: bool = False,
) -> Workflow:
    """
    Submits a conformer-search workflow to the API.

    :param initial_molecule: Molecule to perform the conformer search on.
    :param conf_gen_settings: Conformer generation method and settings. Defaults to
        ``ETKDGSettings()``. Available options (importable directly from ``rowan``):

        - ``ETKDGSettings``  -- RDKit ETKDG, fast, good for most small molecules
        - ``LyrebirdSettings``  -- Rowan ML model
        - ``iMTDGCSettings``  -- CREST iMTD-GC metadynamics, more thorough
        - ``MonteCarloMultipleMinimumSettings``  -- MCMM conformer search
    :param final_method: Method to use for the final optimization.
    :param solvent: Solvent to use for the final optimization.
    :param transition_state: Whether to optimize the transition state.
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
    if conf_gen_settings is None:
        conf_gen_settings = stjames.ETKDGSettings()

    mol_dict = molecule_to_dict(initial_molecule)

    if isinstance(final_method, str):
        final_method = stjames.Method(final_method)

    solvent_model = None
    if solvent:
        solvent_model = "alpb" if final_method in stjames.XTB_METHODS else "cpcm"

    opt_settings = stjames.Settings(
        method=final_method,
        tasks=["optimize"],
        mode=stjames.Mode.AUTO,
        solvent_settings={"solvent": solvent, "model": solvent_model} if solvent else None,
        opt_settings={"transition_state": transition_state, "constraints": []},
    )

    msos = stjames.MultiStageOptSettings(
        mode=stjames.Mode.MANUAL,
        xtb_preopt=True,
        optimization_settings=[opt_settings],
    )

    workflow = stjames.ConformerSearchWorkflow(
        initial_molecule=mol_dict,
        multistage_opt_settings=msos,
        conf_gen_settings=conf_gen_settings,
        solvent=solvent,
        transition_state=transition_state,
    )

    data = {
        "workflow_type": "conformer_search",
        "workflow_data": workflow.model_dump(mode="json"),
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
