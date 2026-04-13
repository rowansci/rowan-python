"""Basic calculation workflow - perform quantum chemical calculations."""

from typing import Any, Literal, TypedDict, cast

import stjames
from stjames import (
    BasisSet,
    Correction,
    Engine,
    Method,
    OptimizationSettings,
    PBCDFTSettings,
    Settings,
    SolventSettings,
)

from ..calculation import Calculation, retrieve_calculation
from ..folder import Folder
from ..molecule import Molecule
from ..utils import api_client
from .base import (
    MoleculeInput,
    Workflow,
    WorkflowResult,
    molecule_to_dict,
    register_result,
)
from .constants import to_relative_kcal

PresetName = Literal[
    "general_nnp", "organic_nnp", "rapid_semiempirical", "routine_dft", "careful_dft"
]


class _PresetsDict(TypedDict, total=True):
    general_nnp: Settings
    organic_nnp: Settings
    rapid_semiempirical: Settings
    routine_dft: Settings
    careful_dft: Settings


# Named presets for common level-of-theory combinations
_PRESETS: _PresetsDict = {
    "general_nnp": Settings(method=Method.OMOL25_CONSERVING_S),
    "organic_nnp": Settings(method=Method.AIMNET2_WB97MD3),
    "rapid_semiempirical": Settings(method=Method.GFN2_XTB),
    "routine_dft": Settings(method=Method.R2SCAN, basis_set="vDZP", corrections=[Correction.D4]),
    "careful_dft": Settings(method=Method.WB97MD3BJ, basis_set="vDZP"),
}


@register_result("basic_calculation")
class BasicCalculationResult(WorkflowResult):
    """Result from a basic-calculation workflow."""

    _stjames_class = stjames.BasicCalculationWorkflow

    def __post_init__(self) -> None:
        super().__post_init__()
        if self.complete:
            if calc_uuid := getattr(self._workflow, "calculation_uuid", None):
                self._cache["calculation"] = retrieve_calculation(calc_uuid)

    def __repr__(self) -> str:
        energy = self.energy
        e_str = f"{energy:.6f}" if energy is not None else "None"
        return f"<BasicCalculationResult energy={e_str} H>"

    @property
    def calculation_uuid(self) -> str | None:
        """The UUID of the calculation."""
        return getattr(self._workflow, "calculation_uuid", None)

    @property
    def calculation(self) -> Calculation | None:
        """Lazily fetched Calculation object with full molecule data."""
        if "calculation" not in self._cache:
            if calc_uuid := self.calculation_uuid:
                self._cache["calculation"] = retrieve_calculation(calc_uuid)
        return self._cache.get("calculation")

    @property
    def molecule(self) -> Molecule | None:
        """The final molecule geometry with all computed properties."""
        calc = self.calculation
        return calc.molecule if calc else None

    @property
    def molecules(self) -> list[Molecule]:
        """All molecules from the calculation (e.g., optimization trajectory)."""
        calc = self.calculation
        return calc.molecules if calc else []

    @property
    def energy(self) -> float | None:
        """Energy of the final molecule (Hartree)."""
        mol = self.molecule
        return mol.energy if mol else None

    def optimization_energies(self, relative: bool = False) -> list[float]:
        """
        Energies for each optimization step.

        :param relative: If True, return relative energies in kcal/mol (relative to
            the lowest energy step). If False (default), return absolute energies
            in Hartree.
        :returns: List of energies for each optimization step.
        """
        energies: list[float] = [m.energy for m in self.molecules if m.energy is not None]
        return to_relative_kcal(energies) if relative else energies

    @property
    def charges(self) -> list[float] | None:
        """Partial charges on each atom."""
        mol = self.molecule
        return mol.charges if mol else None

    @property
    def spin_densities(self) -> list[float] | None:
        """Spin densities on each atom (for open-shell systems)."""
        mol = self.molecule
        return mol.spin_densities if mol else None

    @property
    def dipole(self) -> tuple[float, float, float] | None:
        """Dipole moment vector (Debye)."""
        mol = self.molecule
        return mol.dipole if mol else None

    @property
    def frequencies(self) -> list[float] | None:
        """Vibrational frequencies in cm^-1 (if frequency calculation was performed)."""
        mol = self.molecule
        return mol.frequencies if mol else None


def settings_from_preset(preset: PresetName, **overrides: Any) -> stjames.Settings:
    """
    Construct a `Settings` object from a named preset.

    :param preset: Preset name — see `PresetName` for options.
    :param overrides: Any `Settings` fields to override (e.g. tasks, mode, solvent_settings).
    :returns: Validated `Settings` object.
    :raises ValueError: if an unknown preset name is given.
    """
    if preset not in _PRESETS:
        raise ValueError(f"Unknown preset {preset!r}. Choose from: {list(_PRESETS)}")
    base = _PRESETS[cast(PresetName, preset)]
    return Settings.model_validate({**base.model_dump(), **overrides})


def submit_basic_calculation_workflow(
    initial_molecule: MoleculeInput,
    tasks: list[str],
    method: Method | str | None = None,
    basis_set: BasisSet | str | None = None,
    mode: str = "auto",
    engine: Engine | str | None = None,
    corrections: list[str] | None = None,
    solvent_settings: SolventSettings | dict[str, Any] | None = None,
    opt_settings: OptimizationSettings | dict[str, Any] | None = None,
    pbc_dft_settings: PBCDFTSettings | dict[str, Any] | None = None,
    preset: PresetName | None = None,
    name: str = "Basic Calculation Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
    is_draft: bool = False,
) -> Workflow:
    """
    Submit a basic-calculation workflow to the API.

    :param initial_molecule: Molecule to perform the calculation on.
    :param tasks: Tasks to perform, see `Task` (e.g. optimize, energy, frequencies).
    :param method: Computational method, see `Method`. Default: `omol25_conserving_s`.
    :param basis_set: Basis set, see `BasisSet`.
    :param mode: Accuracy mode, see `Mode`. Default: `auto`.
    :param engine: Compute engine, see `Engine`. Auto-selected from method if not specified.
    :param corrections: Dispersion corrections, see `Correction`.
    :param solvent_settings: Solvent settings as a dict or `SolventSettings`.
    :param opt_settings: Optimization settings as a dict or `OptimizationSettings`.
    :param pbc_dft_settings: Periodic boundary condition DFT settings as a dict or
        `PBCDFTSettings`. Specifies the plane-wave cutoff (Hartree), Monkhorst–Pack
        k-point grid, and optional smearing. When set, the engine is automatically
        set to Quantum ESPRESSO unless ``engine`` is explicitly provided.
    :param preset: Named preset, mutually exclusive with method/engine/basis_set/corrections.
        - `general_nnp` — omol25_conserving_s on omol25
        - `organic_nnp` — aimnet2_wb97md3 on aimnet2
        - `rapid_semiempirical` — gfn2_xtb on xtb
        - `routine_dft` — r2scan-D4/vDZP on gpu4pyscf
        - `careful_dft` — wb97m_d3bj/vDZP on gpu4pyscf
    :param name: Workflow name.
    :param folder_uuid: UUID of the folder to place the workflow in.
    :param folder: Folder to place the workflow in.
    :param max_credits: Maximum credits to use.
    :param webhook_url: URL that Rowan will POST to when the workflow completes.
    :param is_draft: If True, submit the workflow as a draft without starting execution.
    :returns: Submitted workflow.
    :raises requests.HTTPError: if the API request fails.
    :raises ValueError: if preset is combined with method/engine/basis_set/corrections.

    To resubmit with identical settings insulated from future preset changes, use
    `submit_workflow` directly with `workflow_data` from a previous result.
    """
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder:
        folder_uuid = folder.uuid

    if isinstance(pbc_dft_settings, dict):
        pbc_dft_settings = PBCDFTSettings(**pbc_dft_settings)

    if preset is not None:
        if any(x is not None for x in [method, engine, corrections, basis_set]):
            raise ValueError(
                "preset is mutually exclusive with method, engine, basis_set, and corrections."
            )
        settings = settings_from_preset(
            preset,
            tasks=tasks,
            mode=mode,
            solvent_settings=solvent_settings,
            **({"opt_settings": opt_settings} if opt_settings else {}),
            **({"pbc_dft_settings": pbc_dft_settings} if pbc_dft_settings else {}),
        )
    else:
        if method is None:
            method = "omol25_conserving_s"

        if isinstance(method, str):
            method = stjames.Method(method)

        if isinstance(solvent_settings, dict):
            solvent_settings = stjames.SolventSettings(**solvent_settings)

        if isinstance(opt_settings, dict):
            opt_settings = stjames.OptimizationSettings(**opt_settings)

        # pbc_dft_settings implies Quantum ESPRESSO; override engine unless explicitly set
        if pbc_dft_settings is not None and engine is None:
            engine = Engine.QUANTUM_ESPRESSO

        settings_kwargs: dict[str, Any] = {
            "method": method,
            "basis_set": basis_set,
            "tasks": tasks,
            "mode": mode,
            "corrections": corrections or [],
            "solvent_settings": solvent_settings,
        }
        if engine is not None:
            settings_kwargs["engine"] = engine
        if opt_settings is not None:
            settings_kwargs["opt_settings"] = opt_settings
        if pbc_dft_settings is not None:
            settings_kwargs["pbc_dft_settings"] = pbc_dft_settings
        settings = stjames.Settings(**settings_kwargs)

    initial_molecule = molecule_to_dict(initial_molecule)

    workflow = stjames.BasicCalculationWorkflow(
        initial_molecule=initial_molecule,
        settings=settings,
        tasks=settings.tasks,
        engine=settings.engine,
    )

    data = {
        "workflow_type": "basic_calculation",
        "workflow_data": workflow.model_dump(mode="json"),
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
