"""Conformer-search workflow - find low-energy molecular conformations."""

import math
from collections.abc import Sequence
from typing import Any

import stjames

from rowan.calculation import Calculation, retrieve_calculation
from rowan.folder import Folder
from rowan.molecule import Molecule
from rowan.types import SolventInput
from rowan.utils import api_client

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
from .constants import ENSEMBLE_TEMPERATURE, boltzmann_weights, to_relative_kcal

ConformerProperties = stjames.ConformerProperties


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
        """Get conformer energies.

        Args:
            relative: return relative energies in kcal/mol (relative to
                the lowest energy conformer). If False (default), return absolute
                energies in Hartree

        Returns:
            list of conformer energies ordered by energy (lowest first)
        """
        energies: list[float] = list(self._workflow.energies)
        return to_relative_kcal(energies) if relative else energies

    @property
    def radii_of_gyration(self) -> list[float | None]:
        """Radius of gyration for each conformer (A), None where it failed."""
        return [p.radius_of_gyration for p in self._workflow.conformer_properties]

    @property
    def sasa(self) -> list[float | None]:
        """Solvent accessible surface area for each conformer (A^2), None where it failed."""
        return [p.solvent_accessible_surface_area for p in self._workflow.conformer_properties]

    @property
    def polar_sasa(self) -> list[float | None]:
        """Polar solvent accessible surface area for each conformer (A^2), None where it failed."""
        return [
            p.polar_solvent_accessible_surface_area for p in self._workflow.conformer_properties
        ]

    @property
    def dipoles(self) -> list[float | None]:
        """Dipole moment magnitude about the center of mass for each conformer (Debye).

        None where the method reports no dipole.
        """
        return [p.dipole for p in self._workflow.conformer_properties]

    @property
    def ensemble_properties(self) -> ConformerProperties | None:
        """Boltzmann-weighted ensemble properties at 298.15 K, as computed by Rowan.

        See `get_ensemble_properties` for the fields and for other temperatures.
        """
        return self._workflow.ensemble_properties

    def get_boltzmann_weights(self, temperature: float = ENSEMBLE_TEMPERATURE) -> list[float]:
        """Compute Boltzmann weights of the conformers from their energies.

        Args:
            temperature: temperature in Kelvin

        Returns:
            weight of each conformer (sum to 1.0), ordered by energy

        Raises:
            ValueError: temperature is not positive or there are no conformers
        """
        return boltzmann_weights(list(self._workflow.energies), temperature)

    def get_ensemble_properties(
        self, temperature: float = ENSEMBLE_TEMPERATURE
    ) -> ConformerProperties:
        """Boltzmann-weighted properties of the conformer ensemble.

        Each property is averaged over the conformers that have it. The dipole is the
        root-mean-square magnitude, and Mulliken charges are averaged per atom. Spectra are
        stick spectra holding every conformer's lines with weighted intensities:
        `excitation_energies` (Hartree) with `oscillator_strengths` (UV-vis) and
        `rotatory_strengths` (ECD), and `vibrational_frequencies` (cm^-1) with
        `ir_intensities` (km/mol). Pass them to `rowan.broaden_spectrum` for the averaged
        spectrum. Rotatory strengths are omitted when the conformers differ in
        stereochemistry, since mixing enantiomers cancels the ECD signal.

        At 298.15 K this returns `ensemble_properties`, computed by Rowan. At other
        temperatures it reweights the per-conformer values. When the ensemble has charges
        or spectra, those are read from every conformer's calculation, which makes one API
        call per conformer on first access.

        Args:
            temperature: temperature in Kelvin

        Returns:
            ensemble properties

        Raises:
            ValueError: temperature is not positive or there are no conformers
        """
        weights = self.get_boltzmann_weights(temperature)
        reference = self._workflow.ensemble_properties
        if temperature == ENSEMBLE_TEMPERATURE and reference is not None:
            return reference

        conformer_props = self._workflow.conformer_properties
        if len(conformer_props) != len(weights):
            conformer_props = [ConformerProperties()] * len(weights)

        mean_square_dipole = _weighted_mean(
            weights, [None if p.dipole is None else p.dipole**2 for p in conformer_props]
        )
        properties: dict[str, Any] = {
            "solvent_accessible_surface_area": _weighted_mean(
                weights, [p.solvent_accessible_surface_area for p in conformer_props]
            ),
            "polar_solvent_accessible_surface_area": _weighted_mean(
                weights, [p.polar_solvent_accessible_surface_area for p in conformer_props]
            ),
            "radius_of_gyration": _weighted_mean(
                weights, [p.radius_of_gyration for p in conformer_props]
            ),
            "dipole": None if mean_square_dipole is None else math.sqrt(mean_square_dipole),
        }
        if reference is not None and (
            reference.mulliken_charges
            or reference.excitation_energies
            or reference.vibrational_frequencies
        ):
            molecules = self._final_molecules()
            excitations, oscillator, rotatory = _pooled_excitations(
                weights, molecules, include_rotatory=reference.rotatory_strengths is not None
            )
            frequencies, ir_intensities = _pooled_ir(weights, molecules)
            properties |= {
                "mulliken_charges": _weighted_mean_lists(
                    weights, [None if m is None else m.mulliken_charges for m in molecules]
                ),
                "excitation_energies": excitations,
                "oscillator_strengths": oscillator,
                "rotatory_strengths": rotatory,
                "vibrational_frequencies": frequencies,
                "ir_intensities": ir_intensities,
            }
        return ConformerProperties(**properties)

    def _final_molecules(self) -> list[stjames.Molecule | None]:
        """Final molecule of each conformer's calculation, None where unavailable."""
        molecules: list[stjames.Molecule | None] = []
        for i in range(self.num_conformers):
            try:
                molecule = self.get_conformer(i).molecule
            except ValueError:
                molecule = None
            molecules.append(None if molecule is None else molecule.to_stjames())
        return molecules

    def get_conformers(self, n: int | None = None) -> list[Molecule]:
        """Fetch conformer molecules.

        Args:
            n: number of conformers to fetch (default: all). Conformers are
                ordered by energy, so n=5 returns the 5 lowest-energy conformers

        Returns:
            molecules

        Note:
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
        """Fetch a conformer's calculation data by index.

        Note:
            Makes one API call per conformer on first access.
            Results are cached. Call clear_cache() to refresh.

        Args:
            index: conformer index (0-based)
            stage: optimization stage (-1 for final stage)

        Returns:
            calculation object with molecule and energy data

        Raises:
            IndexError: index is out of range
            ValueError: conformer UUID is None
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


def _renormalized(weights: Sequence[float], present: Sequence[bool]) -> list[float] | None:
    """Weights restricted to conformers where `present` holds, renormalized to sum to 1.

    Args:
        weights: Boltzmann weight of every conformer
        present: whether each conformer has the property

    Returns:
        renormalized weights with zeros for absent conformers, or None if none is present

    >>> _renormalized([0.5, 0.25, 0.25], [True, False, True])
    [0.6666666666666666, 0.0, 0.3333333333333333]
    """
    total = sum(w for w, p in zip(weights, present, strict=True) if p)
    if not total:
        return None
    return [w / total if p else 0.0 for w, p in zip(weights, present, strict=True)]


def _weighted_mean(weights: Sequence[float], values: Sequence[float | None]) -> float | None:
    """Weighted mean over the conformers that have a value.

    Args:
        weights: Boltzmann weight of every conformer
        values: value of each conformer, None where missing

    Returns:
        weighted mean, or None if no conformer has a value

    >>> _weighted_mean([0.75, 0.25], [1.0, None])
    1.0
    """
    w = _renormalized(weights, [v is not None for v in values])
    if w is None:
        return None
    return sum(wi * v for wi, v in zip(w, values, strict=True) if v is not None)


def _weighted_mean_lists(
    weights: Sequence[float], values: Sequence[Sequence[float] | None]
) -> list[float] | None:
    """Element-wise weighted mean of per-conformer lists over the conformers that have one.

    Args:
        weights: Boltzmann weight of every conformer
        values: equal-length list of each conformer, None where missing

    Returns:
        weighted mean list, or None if no conformer has one

    >>> [round(q, 6) for q in _weighted_mean_lists([0.75, 0.25], [[0.2, -0.2], [0.6, -0.6]])]
    [0.3, -0.3]
    """
    w = _renormalized(weights, [v is not None for v in values])
    if w is None:
        return None
    present = [(wi, v) for wi, v in zip(w, values, strict=True) if v is not None]
    return [sum(wi * v[i] for wi, v in present) for i in range(len(present[0][1]))]


def _pooled_excitations(
    weights: Sequence[float],
    molecules: Sequence[stjames.Molecule | None],
    include_rotatory: bool,
) -> tuple[list[float] | None, list[float] | None, list[float] | None]:
    """Pool every conformer's excitations with weighted oscillator and rotatory strengths.

    Args:
        weights: Boltzmann weight of every conformer
        molecules: final molecule of each conformer, None where unavailable
        include_rotatory: whether to pool rotatory strengths

    Returns:
        (excitation energies in Hartree, weighted oscillator strengths, weighted rotatory
        strengths in 1e-40 esu^2 cm^2); each None when unavailable. Rotatory strengths
        require every conformer with excited states to have them
    """
    has_states = [
        m is not None
        and m.excited_state_energies is not None
        and m.oscillator_strengths is not None
        for m in molecules
    ]
    w = _renormalized(weights, has_states)
    if w is None:
        return None, None, None

    energies: list[float] = []
    oscillator: list[float] = []
    rotatory: list[float] | None = [] if include_rotatory else None
    for wi, m in zip(w, molecules, strict=True):
        if m is None or m.excited_state_energies is None or m.oscillator_strengths is None:
            continue
        ground, *excited = m.excited_state_energies
        energies += [e - ground for e in excited]
        oscillator += [wi * f for f in m.oscillator_strengths[1:]]
        if rotatory is not None:
            if m.rotatory_strengths is None:
                rotatory = None
            else:
                rotatory += [wi * r for r in m.rotatory_strengths[1:]]
    return energies, oscillator, rotatory


def _pooled_ir(
    weights: Sequence[float], molecules: Sequence[stjames.Molecule | None]
) -> tuple[list[float] | None, list[float] | None]:
    """Pool every conformer's real IR lines with weighted intensities.

    Conformers without IR intensities are excluded and the weights renormalized over the
    rest; imaginary modes are dropped.

    Args:
        weights: Boltzmann weight of every conformer
        molecules: final molecule of each conformer, None where unavailable

    Returns:
        (frequencies in cm^-1, weighted IR intensities in km/mol), or (None, None) when no
        conformer has IR intensities
    """
    has_ir = [
        m is not None
        and bool(m.vibrational_modes)
        and all(mode.ir_intensity is not None for mode in m.vibrational_modes or [])
        for m in molecules
    ]
    w = _renormalized(weights, has_ir)
    if w is None:
        return None, None

    frequencies: list[float] = []
    intensities: list[float] = []
    for wi, m, ok in zip(w, molecules, has_ir, strict=True):
        if not ok or m is None:
            continue
        for mode in m.vibrational_modes or []:
            if mode.frequency > 0 and mode.ir_intensity is not None:
                frequencies.append(mode.frequency)
                intensities.append(wi * mode.ir_intensity)
    return frequencies, intensities


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
    frequencies: bool = False,
) -> MultiStageOptSettings:
    """Build optimization settings for the final method.

    Use conformer-search presets, or a single optimization stage for methods
    without a named preset.

    `solvent` and `transition_state` are merged into each stage's `Settings`
    when supplied, and `frequencies` runs on the optimization stage; the named
    preset itself is independent of them (mirrors how tinbergen's form keeps the
    preset dropdown separate from the solvent / TS toggles).
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
        frequencies=frequencies,
    )


def _release_ts_optimization_constraints(
    multistage_opt_settings: MultiStageOptSettings,
) -> MultiStageOptSettings:
    """Remove constraints from transition-state refinement stages."""
    optimization_settings = [
        stage.model_copy(
            update={
                "opt_settings": stage.opt_settings.model_copy(update={"constraints": ()}),
            }
        )
        if stage.opt_settings.transition_state and stage.opt_settings.constraints
        else stage
        for stage in multistage_opt_settings.optimization_settings
    ]
    return multistage_opt_settings.model_copy(
        update={"optimization_settings": optimization_settings}
    )


def submit_conformer_search_workflow(
    initial_molecule: StructureInput | SMILES | None = None,
    conf_gen_settings: ConformerGenSettings | None = None,
    final_method: Method | str = "aimnet2_wb97md3",
    solvent: SolventInput = None,
    transition_state: bool = False,
    frequencies: bool = False,
    multistage_opt_settings: MultiStageOptSettings | None = None,
    conformer_clustering_settings: ConformerClusteringSettings | None = None,
    initial_conformers: Sequence[StructureInput] | None = None,
    name: str = "Conformer Search Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
    is_draft: bool = False,
) -> Workflow[ConformerSearchResult]:
    """Submits a conformer-search workflow to the API.

    Runs in one of two modes:

    - **Generate** (default): build conformers from `initial_molecule` using
      `conf_gen_settings`, then optimize, deduplicate, and rank them.
    - **Screen-only**: pass `initial_conformers` (and leave `conf_gen_settings` as
      `None`) to skip generation and run only optimize / deduplicate / rank on
      conformers you already have. Useful when geometries come from another tool
      (RDKit, CREST, OMEGA), a crystal or MD ensemble, or a previous workflow, and
      you want consistent optimized energies and a deduplicated ranked ensemble.

    Args:
        initial_molecule: molecule to perform the conformer search on (omit when using
            `initial_conformers`). A 3D structure for any generator; a SMILES string is also
            accepted when `conf_gen_settings` is `ETKDGSettings` or `OpenConfSettings`
            (which build geometry from topology)
        conf_gen_settings: conformer generation method and settings. Defaults to
            `OpenConfSettings()`. Available options (importable directly from `rowan`):

            - `ETKDGSettings`  -- RDKit ETKDG, fast, good for most small molecules (SMILES ok)
            - `OpenConfSettings`  -- OpenConf generator (SMILES ok)
            - `iMTDGCSettings`  -- CREST iMTD-GC metadynamics, more thorough (3D structure only)
            - `MonteCarloMultipleMinimumSettings`  -- MCMM conformer search (3D structure only)
        final_method: method to use for the final optimization. Ignored if
            `multistage_opt_settings` is provided
        solvent: solvent to use for the final optimization. Ignored if
            `multistage_opt_settings` is provided
        transition_state: whether to optimize the transition state. Ignored
            if `multistage_opt_settings` is provided
        frequencies: compute vibrational frequencies of each conformer at the optimization
            level, for the ensemble IR spectrum. Ignored if `multistage_opt_settings` is
            provided; set its `frequencies` instead
        multistage_opt_settings: optimization stages and singlepoint settings
            for ranking conformers. When provided, takes precedence over
            `final_method` / `solvent` / `transition_state`. When omitted, an MSO is
            built from those three params. Constraints remain active during conformer
            generation, but are omitted from transition-state refinement stages
        conformer_clustering_settings: cluster the generated ensemble (ReSCoSS k-means on
            3D-shape descriptors) and keep only representative conformers for downstream
            optimization. Not supported with `initial_conformers`
        initial_conformers: pre-generated 3D conformers to optimize, deduplicate, and
            rank directly, skipping conformer generation (screen-only mode). Requirements
            (all enforced):

            - mutually exclusive with `initial_molecule`
            - `conf_gen_settings` must be `None`
            - every conformer must be a real 3D structure (no SMILES)
            - every conformer must be the same molecule with **identical atom ordering** --
              conformers are compared atom-by-atom during deduplication, so atom *i* must be
              the same atom in every conformer. Read them from one multi-conformer source
              (one RDKit mol, an SDF, an MD trajectory) rather than assembling them separately
        name: name of the workflow
        folder_uuid: UUID of the folder to place the workflow in
        folder: destination folder
        max_credits: maximum credits for the workflow
        webhook_url: URL that Rowan will POST to when the workflow completes
        is_draft: save as a draft without starting execution

    Returns:
        submitted workflow

    Raises:
        httpx.HTTPStatusError: request to the API fails
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
            final_method,
            solvent=solvent,
            transition_state=transition_state,
            frequencies=frequencies,
        )
    multistage_opt_settings = _release_ts_optimization_constraints(multistage_opt_settings)

    workflow = stjames.ConformerSearchWorkflow(
        initial_molecule=mol_dict,
        initial_smiles=initial_smiles,
        initial_conformers=conformer_dicts,
        multistage_opt_settings=multistage_opt_settings,
        conf_gen_settings=conf_gen_settings,
        conformer_clustering_settings=conformer_clustering_settings,
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
        return Workflow[ConformerSearchResult](**response.json())
