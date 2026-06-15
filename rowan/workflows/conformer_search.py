"""Conformer-search workflow - find low-energy molecular conformations."""

from typing import Any

import stjames

from ..calculation import Calculation, retrieve_calculation
from ..folder import Folder
from ..molecule import Molecule
from ..types import SolventInput
from ..utils import api_client
from .base import (
    SMILES,
    ConformerClusteringSettings,
    ConformerGenSettings,
    ETKDGSettings,
    Method,
    Mode,
    MultiStageOptSettings,
    OpenConfSettings,
    Settings,
    StructureInput,
    Task,
    Workflow,
    WorkflowResult,
    molecule_to_dict,
    register_result,
    require_coordinates,
)
from .constants import to_relative_kcal


@register_result("conformer_search")
class ConformerSearchResult(WorkflowResult):
    """Result from a conformer-search workflow."""

    _stjames_class = stjames.ConformerSearchWorkflow

    def __repr__(self) -> str:
        n = self.num_conformers
        return f"<ConformerSearchResult conformers={n}>"

    def __post_init__(self) -> None:
        """Default `conf_gen_settings` to None (omitted by the server in skip mode), then parse."""
        self.workflow_data.setdefault("conf_gen_settings", None)
        super().__post_init__()

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


_FINAL_METHOD_TO_MSO_RECIPE: dict[Method, tuple[Method, Method | None]] = {
    Method.AIMNET2_WB97MD3: (Method.AIMNET2_WB97MD3, None),
    Method.OMOL25_CONSERVING_S: (Method.OMOL25_CONSERVING_S, None),
    Method.GFN2_XTB: (Method.GFN2_XTB, None),
    Method.G_XTB: (Method.GFN2_XTB, Method.G_XTB),
    Method.R2SCAN3C: (Method.GFN2_XTB, Method.R2SCAN3C),
}


def _mso_for_final_method(
    final_method: Method,
    solvent: SolventInput = None,
    transition_state: bool = False,
) -> MultiStageOptSettings:
    """Build a `MultiStageOptSettings` for `final_method` matching tinbergen's
    conformer-search MSO presets. Falls back to a single opt stage at
    `final_method` for methods without a named preset.

    `solvent` and `transition_state` are merged into each stage's `Settings`
    when supplied; the named preset itself is independent of them (mirrors how
    tinbergen's form keeps the preset dropdown separate from the solvent / TS
    toggles).
    """
    opt_method, sp_method = _FINAL_METHOD_TO_MSO_RECIPE.get(final_method, (final_method, None))

    def _solvent_settings_for(method: Method) -> dict | None:
        if not solvent:
            return None
        return {
            "solvent": solvent,
            "model": "alpb" if method in stjames.XTB_METHODS else "cpcm",
        }

    opt_stage = Settings(
        method=opt_method,
        tasks=[Task.OPTIMIZE],
        mode=Mode.AUTO,
        solvent_settings=_solvent_settings_for(opt_method),
        opt_settings={"transition_state": transition_state, "constraints": []},
    )
    sp_stage: Settings | None = None
    if sp_method is not None:
        sp_stage = Settings(
            method=sp_method,
            tasks=[Task.ENERGY],
            mode=Mode.AUTO,
            solvent_settings=_solvent_settings_for(sp_method),
        )

    return MultiStageOptSettings(
        optimization_settings=[opt_stage],
        singlepoint_settings=sp_stage,
    )


def submit_conformer_search_workflow(
    initial_molecule: StructureInput | SMILES | None = None,
    conf_gen_settings: ConformerGenSettings | None = None,
    final_method: Method | str = "aimnet2_wb97md3",
    solvent: SolventInput = None,
    transition_state: bool = False,
    multistage_opt_settings: MultiStageOptSettings | None = None,
    conformer_clustering_settings: ConformerClusteringSettings | None = None,
    initial_conformers: list[StructureInput] | None = None,
    name: str = "Conformer Search Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
    is_draft: bool = False,
) -> Workflow:
    """
    Submits a conformer-search workflow to the API.

    Runs in one of two modes:

    - **Generate** (default): build conformers from `initial_molecule` using
      `conf_gen_settings`, then optimize, deduplicate, and rank them.
    - **Screen-only**: pass `initial_conformers` (and leave `conf_gen_settings` as
      ``None``) to skip generation and run only optimize / deduplicate / rank on
      conformers you already have. Useful when geometries come from another tool
      (RDKit, CREST, OMEGA), a crystal or MD ensemble, or a previous workflow, and
      you want consistent optimized energies and a deduplicated ranked ensemble.

    :param initial_molecule: Molecule to perform the conformer search on (omit when using
        `initial_conformers`). A 3D structure for any generator; a SMILES string is also
        accepted when `conf_gen_settings` is ``ETKDGSettings`` or ``OpenConfSettings``
        (which build geometry from topology).
    :param conf_gen_settings: Conformer generation method and settings. Defaults to
        ``OpenConfSettings()``. Available options (importable directly from ``rowan``):

        - ``ETKDGSettings``  -- RDKit ETKDG, fast, good for most small molecules (SMILES ok)
        - ``OpenConfSettings``  -- OpenConf generator (SMILES ok)
        - ``iMTDGCSettings``  -- CREST iMTD-GC metadynamics, more thorough (3D structure only)
        - ``MonteCarloMultipleMinimumSettings``  -- MCMM conformer search (3D structure only)
    :param final_method: Method to use for the final optimization. Ignored if
        `multistage_opt_settings` is provided.
    :param solvent: Solvent to use for the final optimization. Ignored if
        `multistage_opt_settings` is provided.
    :param transition_state: Whether to optimize the transition state. Ignored
        if `multistage_opt_settings` is provided.
    :param multistage_opt_settings: Optimization stages and singlepoint settings
        for ranking conformers. When provided, takes precedence over
        `final_method` / `solvent` / `transition_state`. When omitted, an MSO is
        built from those three params.
    :param conformer_clustering_settings: Cluster the generated ensemble (ReSCoSS k-means on
        3D-shape descriptors) and keep only representative conformers for downstream
        optimization. Not supported with `initial_conformers`.
    :param initial_conformers: Pre-generated 3D conformers to optimize, deduplicate, and
        rank directly, skipping conformer generation (screen-only mode). Requirements
        (all enforced):

        - mutually exclusive with `initial_molecule`
        - `conf_gen_settings` must be ``None``
        - every conformer must be a real 3D structure (no SMILES)
        - every conformer must be the same molecule with **identical atom ordering** --
          conformers are compared atom-by-atom during deduplication, so atom *i* must be
          the same atom in every conformer. Read them from one multi-conformer source
          (one RDKit mol, an SDF, an MD trajectory) rather than assembling them separately.
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

    initial_smiles = ""
    mol_dict: dict[str, Any] | None = None
    conformer_dicts: list[dict[str, Any]] = []

    if initial_conformers:
        # Screen-only mode: rank the supplied conformers directly, no generation.
        if initial_molecule is not None:
            raise ValueError("Provide either `initial_molecule` or `initial_conformers`, not both.")
        if conf_gen_settings is not None:
            raise ValueError(
                "`conf_gen_settings` must be None when using `initial_conformers`; "
                "screen-only mode does not generate conformers."
            )
        if conformer_clustering_settings is not None:
            raise ValueError(
                "`conformer_clustering_settings` is not supported with `initial_conformers`; "
                "clustering selects conformers from a generated ensemble, which screen-only "
                "mode does not produce."
            )
        reference_atoms: list[int] | None = None
        for conformer in initial_conformers:
            require_coordinates(conformer)
            conformer_dict = molecule_to_dict(conformer)
            atoms = [atom["atomic_number"] for atom in conformer_dict["atoms"]]
            if reference_atoms is None:
                reference_atoms = atoms
            elif atoms != reference_atoms:
                raise ValueError(
                    "All `initial_conformers` must be the same molecule with identical atom "
                    "ordering; conformers are compared atom-by-atom, so every conformer must "
                    "list the same atoms in the same order."
                )
            conformer_dicts.append(conformer_dict)
    else:
        if initial_molecule is None:
            raise ValueError("Provide `initial_molecule` or `initial_conformers`.")
        if conf_gen_settings is None:
            conf_gen_settings = OpenConfSettings()
        # ETKDG and OpenConf build geometry from a bare SMILES; the other generators
        # (CREST, MCMM) need a real 3D structure.
        generator_builds_geometry = isinstance(conf_gen_settings, (ETKDGSettings, OpenConfSettings))
        if isinstance(initial_molecule, str):
            if not generator_builds_geometry:
                raise ValueError(
                    "SMILES input is only supported with ETKDG or OpenConf conformer generation. "
                    "Provide a 3D structure (rowan.Molecule.from_smiles(...) or from_xyz(...)) for "
                    f"{type(conf_gen_settings).__name__}."
                )
            initial_smiles = initial_molecule
        else:
            if not generator_builds_geometry:
                require_coordinates(initial_molecule)
            mol_dict = molecule_to_dict(initial_molecule)

    if multistage_opt_settings is None:
        if isinstance(final_method, str):
            final_method = Method(final_method)
        multistage_opt_settings = _mso_for_final_method(
            final_method, solvent=solvent, transition_state=transition_state
        )

    workflow = stjames.ConformerSearchWorkflow(
        initial_molecule=mol_dict,
        initial_smiles=initial_smiles,
        initial_conformers=conformer_dicts,
        multistage_opt_settings=multistage_opt_settings,
        conf_gen_settings=conf_gen_settings,
        conformer_clustering_settings=conformer_clustering_settings,
        solvent=solvent,
        transition_state=transition_state,
    )

    workflow_data = workflow.model_dump(serialize_as_any=True, mode="json")
    # model_dump drops `conf_gen_settings` when it is None (skip-generation mode), but the
    # server requires the key to be present; re-add it as null.
    workflow_data.setdefault("conf_gen_settings", None)
    data = {
        "workflow_type": "conformer_search",
        "workflow_data": workflow_data,
        "name": name,
        "folder_uuid": folder_uuid,
        "max_credits": max_credits,
        "webhook_url": webhook_url,
        "is_draft": is_draft,
    }
    if mol_dict is not None:
        data["initial_molecule"] = mol_dict
    elif conformer_dicts:
        data["initial_molecule"] = conformer_dicts[0]
    else:
        data["initial_smiles"] = initial_smiles

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())
