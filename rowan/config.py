"""Configuration for validating method/engine/basis set/workflow compatibility.

This module provides capability information for blocking invalid combinations
of engines, methods, basis sets, elements, and workflow-specific settings.
"""

from dataclasses import dataclass, field
from typing import Callable, Literal

import stjames

Task = stjames.Task


@dataclass(frozen=True, slots=True)
class MethodSettings:
    """Settings and constraints for a computational method."""

    atoms_supported: set[int] | None = None
    disabled_tasks: set[Task] = field(default_factory=set)
    disable_basis_sets: bool = False
    disable_corrections: bool = False
    disable_optimize_cell: bool = False


_TO_RN = set(range(1, 87))  # H to Rn (1-86)
_TO_LR = set(range(1, 104))  # H to Lr (1-103)

METHOD_SETTINGS: dict[str, MethodSettings] = {
    # Neural Network Potentials
    "aimnet2_wb97md3": MethodSettings(
        atoms_supported={1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 33, 34, 35, 53},
    ),
    "egret_1": MethodSettings(
        atoms_supported={1, 6, 7, 8, 9, 15, 16, 17, 35, 53},
    ),
    "egret_1e": MethodSettings(
        atoms_supported={1, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53},
    ),
    "egret_1t": MethodSettings(
        atoms_supported={1, 6, 7, 8, 9, 15, 16, 17, 35, 53},
    ),
    "omol25_conserving_s": MethodSettings(
        disable_optimize_cell=True,
    ),
    "orb_v3_conservative_inf_omat": MethodSettings(
        disable_corrections=True,
    ),
    # xTB methods
    "gfn_ff": MethodSettings(atoms_supported=_TO_LR),
    "gfn0_xtb": MethodSettings(atoms_supported=_TO_RN),
    "gfn1_xtb": MethodSettings(atoms_supported=_TO_RN),
    "gfn2_xtb": MethodSettings(atoms_supported=_TO_RN),
    "g_xtb": MethodSettings(
        atoms_supported=_TO_LR,
        disabled_tasks={
            Task.CHARGE,
            Task.SPIN_DENSITY,
            Task.DIPOLE,
            Task.HESSIAN,
            Task.FREQUENCIES,
        },
    ),
    # Composite methods
    "hf_3c": MethodSettings(
        atoms_supported=_TO_RN,
        disable_basis_sets=True,
        disable_corrections=True,
    ),
    "b97_3c": MethodSettings(
        atoms_supported=_TO_RN,
        disable_basis_sets=True,
        disable_corrections=True,
    ),
    "b97_d3bj": MethodSettings(disable_corrections=True),
    "r2scan_3c": MethodSettings(
        atoms_supported=_TO_RN,
        disable_basis_sets=True,
        disable_corrections=True,
    ),
    "wb97x_3c": MethodSettings(
        atoms_supported=_TO_RN,
        disable_basis_sets=True,
        disable_corrections=True,
    ),
    # DFT with built-in dispersion
    "wb97x_d3": MethodSettings(disable_corrections=True),
    "wb97m_d3bj": MethodSettings(disable_corrections=True),
    "wb97x_v": MethodSettings(disable_corrections=True),
    "wb97m_v": MethodSettings(disable_corrections=True),
    "dsd_blyp_d3bj": MethodSettings(disable_corrections=True),
    # ML functionals
    "skala": MethodSettings(
        disable_corrections=True,
        disabled_tasks={
            Task.OPTIMIZE,
            Task.OPTIMIZE_TS,
            Task.FREQUENCIES,
            Task.HESSIAN,
            Task.GRADIENT,
        },
    ),
}


@dataclass(frozen=True, slots=True)
class EngineSettings:
    """Settings and constraints for a computational engine."""

    methods: set[str]
    atoms_supported: set[int] | None = None
    disabled_tasks: Callable[[bool], set[Task]] | None = None
    disable_basis_sets: bool = False
    disable_corrections: bool = False
    corrections: set[str] = field(default_factory=set)
    disable_solvents: bool = False
    solvent_models: set[str] = field(default_factory=set)
    disable_open_shell: bool = False
    disable_charge: bool = False
    periodic: bool | Literal["only", "warn"] = False


_MP_TRJ = {i for i in range(1, 96) if i not in {84, 85, 86, 87, 88}}
_OMOL25 = set(range(1, 84))  # H to Bi (1-83)

ENGINE_SETTINGS: dict[str, EngineSettings] = {
    "aimnet2": EngineSettings(
        methods={"aimnet2_wb97md3"},
        atoms_supported={1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 33, 34, 35, 53},
        disabled_tasks=lambda periodic: (
            {Task.DIPOLE, Task.SPIN_DENSITY, Task.FREQUENCIES}
            if periodic
            else {Task.DIPOLE, Task.SPIN_DENSITY}
        ),
        disable_basis_sets=True,
        disable_corrections=True,
        solvent_models={"alpb", "cpcmx"},
        disable_open_shell=True,
        periodic=True,
    ),
    "egret": EngineSettings(
        methods={"egret_1", "egret_1e", "egret_1t"},
        disabled_tasks=lambda periodic: (
            {Task.CHARGE, Task.DIPOLE, Task.SPIN_DENSITY, Task.FREQUENCIES}
            if periodic
            else {Task.CHARGE, Task.DIPOLE, Task.SPIN_DENSITY}
        ),
        disable_basis_sets=True,
        disable_corrections=True,
        solvent_models={"alpb", "cpcmx"},
        disable_open_shell=True,
        disable_charge=True,
        periodic=True,
    ),
    "mace": EngineSettings(
        methods={"mace_mp_0b2_l"},
        atoms_supported=_MP_TRJ,
        disabled_tasks=lambda periodic: (
            {Task.CHARGE, Task.DIPOLE, Task.SPIN_DENSITY, Task.FREQUENCIES}
            if periodic
            else {Task.CHARGE, Task.DIPOLE, Task.SPIN_DENSITY}
        ),
        disable_basis_sets=True,
        corrections={"d3bj"},
        solvent_models={"alpb", "cpcmx"},
        disable_open_shell=True,
        disable_charge=True,
        periodic="warn",
    ),
    "omol25": EngineSettings(
        methods={
            "omol25_conserving_s",
            "uma_s_omol",
            "uma_m_omol",
            "uma_s_omat",
            "uma_m_omat",
            "uma_s_omc",
            "uma_m_omc",
        },
        atoms_supported=_OMOL25,
        disabled_tasks=lambda periodic: (
            {Task.CHARGE, Task.DIPOLE, Task.SPIN_DENSITY, Task.FREQUENCIES}
            if periodic
            else {Task.CHARGE, Task.DIPOLE, Task.SPIN_DENSITY}
        ),
        disable_basis_sets=True,
        disable_corrections=True,
        solvent_models={"alpb", "cpcmx"},
        periodic=True,
    ),
    "orb": EngineSettings(
        methods={"orb_v3_conservative_inf_omat"},
        atoms_supported=_MP_TRJ,
        disabled_tasks=lambda periodic: (
            {Task.CHARGE, Task.DIPOLE, Task.SPIN_DENSITY, Task.FREQUENCIES}
            if periodic
            else {Task.CHARGE, Task.DIPOLE, Task.SPIN_DENSITY}
        ),
        disable_basis_sets=True,
        corrections={"d3"},
        solvent_models={"alpb", "cpcmx"},
        disable_open_shell=True,
        disable_charge=True,
        periodic="warn",
    ),
    "xtb": EngineSettings(
        methods={"gfn2_xtb", "g_xtb", "gfn_ff", "gfn0_xtb", "gfn1_xtb"},
        disabled_tasks=lambda periodic: {Task.SPIN_DENSITY},
        disable_basis_sets=True,
        disable_corrections=True,
        solvent_models={"alpb", "cpcmx"},
    ),
    "tblite": EngineSettings(
        methods={"gfn2_xtb"},
        atoms_supported=_TO_RN,
        disabled_tasks=lambda periodic: {Task.SPIN_DENSITY, Task.FREQUENCIES, Task.OPTIMIZE_TS},
        disable_basis_sets=True,
        disable_corrections=True,
        disable_solvents=True,
        periodic="only",
    ),
    "psi4": EngineSettings(
        methods={
            "hf",
            "pbe",
            "r2scan",
            "tpss",
            "m06l",
            "pbe0",
            "b3lyp",
            "tpssh",
            "m06",
            "m062x",
            "camb3lyp",
            "wb97x_v",
            "wb97x_d3",
            "wb97m_v",
            "wb97m_d3bj",
            "hf_3c",
            "b97_3c",
            "b97_d3bj",
            "r2scan_3c",
            "wb97x_3c",
            "dsd_blyp_d3bj",
        },
        corrections={"d3bj"},
        solvent_models={"cpcm", "iefpcm"},
    ),
    "pyscf": EngineSettings(
        methods={
            "pbe",
            "r2scan",
            "tpss",
            "m06l",
            "pbe0",
            "b3lyp",
            "tpssh",
            "m06",
            "m062x",
            "camb3lyp",
            "wb97x_v",
            "wb97m_v",
            "wb97m_d3bj",
            "hf",
            "skala",
        },
        corrections={"d3bj", "d4"},
        solvent_models={"cpcm", "cosmo", "iefpcm", "smd"},
    ),
    "gpu4pyscf": EngineSettings(
        methods={
            "pbe",
            "r2scan",
            "tpss",
            "m06l",
            "pbe0",
            "b3lyp",
            "tpssh",
            "m06",
            "m062x",
            "camb3lyp",
            "wb97x_v",
            "wb97m_v",
            "wb97m_d3bj",
            "hf",
            "skala",
        },
        corrections={"d3bj", "d4"},
        solvent_models={"cpcm", "cosmo", "iefpcm", "smd"},
    ),
}


_TO_AR = set(range(1, 19))  # H to Ar (1-18)
_TO_KR = set(range(1, 37))  # H to Kr (1-36)
_TO_XE = set(range(1, 55))  # H to Xe (1-54)
_TO_KR_MINUS_K = _TO_KR - {19}
_TO_KR_MINUS_K_CA = _TO_KR - {19, 20}
_SC_TO_ZN = set(range(21, 31))
_TO_KR_MINUS_SC_TO_ZN = _TO_KR - _SC_TO_ZN

BASIS_SET_ATOMS: dict[str, set[int]] = {
    "STO-3G": _TO_XE,
    "STO-6G": _TO_XE,
    "MIDIX": {1, 3, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53},
    "MINIX": _TO_RN,
    "mTZVP": _TO_RN,
    "vDZP": _TO_RN,
    "3-21G": _TO_XE,
    "6-31G": _TO_KR,
    "6-31G(d)": _TO_KR,
    "6-31+G(d)": _TO_AR,
    "6-31++G(d)": _TO_AR,
    "cc-pVDZ": _TO_KR_MINUS_K,
    "cc-pVTZ": _TO_KR_MINUS_K,
    "cc-pVQZ": _TO_KR_MINUS_K,
    "cc-pV5Z": _TO_KR_MINUS_K,
    "aug-cc-pVDZ": _TO_KR_MINUS_K_CA,
    "aug-cc-pVTZ": _TO_KR_MINUS_K_CA,
    "aug-cc-pVQZ": _TO_KR_MINUS_K_CA,
    "aug-cc-pV5Z": _TO_KR_MINUS_K_CA,
    "pc-0": _TO_KR_MINUS_SC_TO_ZN,
    "pc-1": _TO_KR,
    "pc-2": _TO_KR,
    "pc-3": _TO_KR,
    "pc-4": _TO_KR,
    "aug-pc-0": _TO_KR_MINUS_SC_TO_ZN,
    "aug-pc-1": _TO_KR,
    "aug-pc-2": _TO_KR,
    "aug-pc-3": _TO_KR,
    "aug-pc-4": _TO_KR,
    "pcseg-0": _TO_KR,
    "pcseg-1": _TO_KR,
    "pcseg-2": _TO_KR,
    "pcseg-3": _TO_KR,
    "pcseg-4": _TO_KR,
    "aug-pcseg-0": _TO_KR,
    "aug-pcseg-1": _TO_KR,
    "aug-pcseg-2": _TO_KR,
    "aug-pcseg-3": _TO_KR,
    "aug-pcseg-4": _TO_KR,
    "def2-SVP": _TO_RN,
    "def2-SVPD": _TO_RN,
    "def2-TZVP": _TO_RN,
    "def2-TZVPD": _TO_RN,
    "def2-TZVPP": _TO_RN,
    "def2-TZVPPD": _TO_RN,
    "def2-QZVP": _TO_RN,
    "def2-QZVPD": _TO_RN,
    "def2-QZVPP": _TO_RN,
    "def2-QZVPPD": _TO_RN,
}


@dataclass(frozen=True, slots=True)
class ConformerGeneratorSettings:
    """
    Constraints for a conformer generator.

    :param disable_constraints: Geometry constraints not supported.
    :param disable_open_shell: Only closed-shell molecules supported.
    :param disable_ts: Transition state search not supported.
    :param atoms_supported: Supported atomic numbers (None = all).
    :param allowed_engines: Engines allowed for energy evaluation (None = all).
    :param solvent_warning: Warn if generator solvent differs from final opt solvent.
    """

    disable_constraints: bool = False
    disable_open_shell: bool = False
    disable_ts: bool = False
    atoms_supported: set[int] | None = None
    allowed_engines: list[str] | None = None
    solvent_warning: bool = False


# Main group elements for ETKDG/Lyrebird
_MAIN_GROUP = {1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53}

CONFORMER_GENERATOR_SETTINGS: dict[str, ConformerGeneratorSettings] = {
    "etkdg": ConformerGeneratorSettings(
        disable_constraints=True,
        disable_open_shell=True,
        disable_ts=True,
        atoms_supported=_MAIN_GROUP,
    ),
    "imtd": ConformerGeneratorSettings(
        disable_constraints=True,
        disable_open_shell=True,
        disable_ts=True,
        solvent_warning=True,  # Warns if iMTD solvent != final optimization solvent
    ),
    "monte_carlo_multiple_minimum": ConformerGeneratorSettings(
        disable_constraints=True,
        disable_open_shell=True,
        disable_ts=True,
        allowed_engines=["aimnet2", "xtb"],  # MCMM only supports these for energy
    ),
    "lyrebird": ConformerGeneratorSettings(
        disable_constraints=True,
        disable_open_shell=True,
        disable_ts=True,
        atoms_supported=_MAIN_GROUP,
    ),
}


@dataclass(frozen=True, slots=True)
class MDSettings:
    """
    Molecular dynamics parameter constraints.

    :param TIMESTEP_MIN: Minimum timestep in femtoseconds.
    :param TIMESTEP_MAX: Maximum timestep in femtoseconds.
    :param MAX_SAVED_FRAMES: Maximum number of frames that can be saved.
    :param NPT_REQUIRES_PERIODIC: NPT ensemble only valid for periodic systems.
    :param QUASICLASSICAL_REQUIRES_FREQUENCIES: Quasiclassical init needs frequencies.
    :param READ_REQUIRES_VELOCITIES: Read init needs velocities from previous job.
    :param CONFINING_BLOCKED_FOR_PERIODIC: Confining constraints invalid for periodic.
    """

    TIMESTEP_MIN: float = 0.5  # fs
    TIMESTEP_MAX: float = 5.0  # fs
    MAX_SAVED_FRAMES: int = 10000

    # Ensemble constraints
    NPT_REQUIRES_PERIODIC: bool = True

    # Initialization constraints
    QUASICLASSICAL_REQUIRES_FREQUENCIES: bool = True
    READ_REQUIRES_VELOCITIES: bool = True

    # Confining constraint
    CONFINING_BLOCKED_FOR_PERIODIC: bool = True


MD_SETTINGS = MDSettings()


@dataclass(frozen=True, slots=True)
class SolubilityMethodSettings:
    """
    Constraints for solubility prediction methods.

    :param allowed_solvents: Solvents supported (None = all). Use SMILES format.
    :param single_temperature_k: If set, only this temperature allowed.
    :param temperature_range: If set, (min_k, max_k) range allowed.
    """

    allowed_solvents: list[str] | None = None
    single_temperature_k: float | None = None
    temperature_range: tuple[float, float] | None = None


SOLUBILITY_METHOD_SETTINGS: dict[str, SolubilityMethodSettings] = {
    "kingfisher": SolubilityMethodSettings(
        allowed_solvents=["O"],  # Water only (SMILES)
        single_temperature_k=298.15,
    ),
    "esol": SolubilityMethodSettings(
        allowed_solvents=["O"],  # Water only
        single_temperature_k=298.15,
    ),
    "fastsolv": SolubilityMethodSettings(
        allowed_solvents=None,  # Any solvent
        temperature_range=(273.15, 373.15),
    ),
}


@dataclass(frozen=True, slots=True)
class MacroPKaSettings:
    """
    MacroPKa workflow parameter constraints.

    :param PH_WARNING_MIN: pH values below this trigger accuracy warning.
    :param PH_WARNING_MAX: pH values above this trigger accuracy warning.
    :param CHARGE_MIN: Minimum allowed charge state.
    :param CHARGE_MAX: Maximum allowed charge state.
    """

    PH_WARNING_MIN: float = 0.0
    PH_WARNING_MAX: float = 20.0
    CHARGE_MIN: int = -7
    CHARGE_MAX: int = 7


MACROPKA_SETTINGS = MacroPKaSettings()


# Tautomers workflow only supports these elements (no transition metals)
TAUTOMER_ATOMS_SUPPORTED = [1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 33, 34, 35, 53]


@dataclass(frozen=True, slots=True)
class CofoldingModelSettings:
    """
    Constraints for protein cofolding models.

    :param max_sequence_length: Maximum total sequence length.
    :param supports_cyclic_peptides: Whether cyclic peptides are supported.
    :param supports_potentials: Whether potentials can be used.
    :param supports_ligand_binding_affinity: Whether binding affinity available.
    :param ligand_max_atoms: Max atoms per ligand for binding affinity.
    :param ligand_atoms_supported: Supported atoms for ligand binding.
    """

    max_sequence_length: int = 3000
    supports_cyclic_peptides: bool = False
    supports_potentials: bool = True
    supports_ligand_binding_affinity: bool = False
    ligand_max_atoms: int | None = None
    ligand_atoms_supported: set[int] | None = None


# Atoms supported for Boltz-2 ligand binding affinity
_BOLTZ2_LIGAND_ATOMS = {
    *range(1, 21),  # H-Ca
    30,
    33,
    34,
    35,
    36,
    37,
    38,
    47,
    52,
    53,
    54,
    55,
    56,
    85,
    86,
    87,
    88,
}

COFOLDING_MODEL_SETTINGS: dict[str, CofoldingModelSettings] = {
    "boltz-1": CofoldingModelSettings(
        max_sequence_length=3000,
    ),
    "boltz-2": CofoldingModelSettings(
        max_sequence_length=3000,
        supports_cyclic_peptides=True,
        supports_ligand_binding_affinity=True,
        ligand_max_atoms=128,
        ligand_atoms_supported=_BOLTZ2_LIGAND_ATOMS,
    ),
    "chai-1r": CofoldingModelSettings(
        max_sequence_length=2048,
        supports_potentials=False,
    ),
}


@dataclass(frozen=True, slots=True)
class BinderDesignSettings:
    """
    Protein binder design parameter constraints.

    :param MAX_SEQUENCE_LENGTH: Maximum total amino acids.
    :param NUM_DESIGNS_MAX_FREE: Max designs for free tier.
    :param NUM_DESIGNS_MAX_SUBSCRIBED: Max designs for subscribed users.
    :param BUDGET_MUST_BE_LEQ_NUM_DESIGNS: budget cannot exceed num_designs.
    :param REQUIRES_DESIGNABLE_REGIONS: At least one sequence must have digits.
    :param REQUIRES_UNIQUE_CHAIN_IDS: All chain IDs must be unique.
    """

    MAX_SEQUENCE_LENGTH: int = 1500
    NUM_DESIGNS_MAX_FREE: int = 100
    NUM_DESIGNS_MAX_SUBSCRIBED: int = 1000
    BUDGET_MUST_BE_LEQ_NUM_DESIGNS: bool = True
    REQUIRES_DESIGNABLE_REGIONS: bool = True
    REQUIRES_UNIQUE_CHAIN_IDS: bool = True


BINDER_DESIGN_SETTINGS = BinderDesignSettings()


# Vinardo scoring is only supported with AutoDock Vina
VINARDO_SUPPORTED_EXECUTABLES = ["autodock_vina"]

# If conformer search is enabled, optimization must also be enabled
DOCKING_CSEARCH_REQUIRES_OPTIMIZATION = True


def get_engine_settings(engine: str) -> EngineSettings | None:
    """Look up settings for a compute engine.

    :param engine: Engine identifier string.
    :returns: Settings for the engine, or None if not found.
    """
    return ENGINE_SETTINGS.get(engine)


def get_method_settings(method: str) -> MethodSettings | None:
    """Look up settings for a calculation method.

    :param method: Method identifier string.
    :returns: Settings for the method, or None if not found.
    """
    return METHOD_SETTINGS.get(method)


def get_engines_for_method(method: str) -> list[str]:
    """Find all engines that support a given method.

    :param method: Method identifier string.
    :returns: List of engine identifier strings that support the method.
    """
    return [engine for engine, settings in ENGINE_SETTINGS.items() if method in settings.methods]


def get_combined_settings(
    engine: str, method: str, periodic: bool = False
) -> MethodSettings | None:
    """Get combined settings for an engine/method combination.

    Merges engine and method constraints to give the full picture of what's supported.

    :param engine: Engine identifier string.
    :param method: Method identifier string.
    :param periodic: Whether periodic boundary conditions are used.
    :returns: Merged MethodSettings, or None if the combination is not supported.
    """
    engine_settings = ENGINE_SETTINGS.get(engine)
    if engine_settings is None or method not in engine_settings.methods:
        return None

    method_settings = METHOD_SETTINGS.get(method, MethodSettings())

    engine_disabled = (
        engine_settings.disabled_tasks(periodic) if engine_settings.disabled_tasks else set()
    )
    method_disabled = method_settings.disabled_tasks
    combined_disabled = engine_disabled | method_disabled

    atoms: set[int] | None
    if engine_settings.atoms_supported and method_settings.atoms_supported:
        atoms = engine_settings.atoms_supported & method_settings.atoms_supported
    else:
        atoms = engine_settings.atoms_supported or method_settings.atoms_supported

    return MethodSettings(
        atoms_supported=atoms,
        disabled_tasks=combined_disabled,
        disable_basis_sets=(
            engine_settings.disable_basis_sets or method_settings.disable_basis_sets
        ),
        disable_corrections=(
            engine_settings.disable_corrections or method_settings.disable_corrections
        ),
        disable_optimize_cell=method_settings.disable_optimize_cell,
    )


def get_conformer_generator_settings(generator: str) -> ConformerGeneratorSettings | None:
    """Look up settings for a conformer generator.

    :param generator: Conformer generator identifier string.
    :returns: Settings for the generator, or None if not found.
    """
    return CONFORMER_GENERATOR_SETTINGS.get(generator)


def get_solubility_method_settings(method: str) -> SolubilityMethodSettings | None:
    """Look up settings for a solubility prediction method.

    :param method: Solubility method identifier string.
    :returns: Settings for the method, or None if not found.
    """
    return SOLUBILITY_METHOD_SETTINGS.get(method)


def get_cofolding_model_settings(model: str) -> CofoldingModelSettings | None:
    """Look up settings for a protein cofolding model.

    :param model: Cofolding model identifier string.
    :returns: Settings for the model, or None if not found.
    """
    return COFOLDING_MODEL_SETTINGS.get(model)
