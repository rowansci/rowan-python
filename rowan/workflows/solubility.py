"""Solubility workflow - predict molecular solubility in various solvents."""

from dataclasses import dataclass
from typing import Literal

import stjames

from ..utils import api_client
from .base import MoleculeInput, Workflow, WorkflowResult, extract_smiles, register_result

# Common solvents with human-readable names (from tinbergen)
# Users can use these names or provide arbitrary SMILES for fastsolv
COMMON_SOLVENTS: dict[str, str] = {
    "acetone": "CC(C)=O",
    "acetonitrile": "CC#N",
    "benzene": "c1ccccc1",
    "chlorobenzene": "Clc1ccccc1",
    "chloroform": "ClC(Cl)Cl",
    "cyclohexane": "C1CCCCC1",
    "dichloromethane": "ClCCl",
    "1,2-dichloroethane": "ClCCCl",
    "diethyl ether": "CCOCC",
    "dimethylformamide": "CN(C)C=O",
    "dmf": "CN(C)C=O",
    "dimethylacetamide": "CN(C)C(C)=O",
    "dma": "CN(C)C(C)=O",
    "dimethylsulfoxide": "CS(C)=O",
    "dmso": "CS(C)=O",
    "1,4-dioxane": "C1COCCO1",
    "dioxane": "C1COCCO1",
    "ethanol": "CCO",
    "ethyl acetate": "CC(=O)OCC",
    "heptane": "CCCCCCC",
    "hexane": "CCCCCC",
    "isopropanol": "CC(C)O",
    "isopropyl acetate": "CC(=O)OC(C)C",
    "methanol": "CO",
    "n-methylpyrrolidone": "CN1CCCC1=O",
    "nmp": "CN1CCCC1=O",
    "nonane": "CCCCCCCCC",
    "pentanol": "CCCCCO",
    "1-pentanol": "CCCCCO",
    "tert-butanol": "CC(C)(C)O",
    "tetrahydrofuran": "C1CCCO1",
    "thf": "C1CCCO1",
    "tetrachloroethylene": "ClC(Cl)=C(Cl)Cl",
    "trichloroethylene": "ClC=C(Cl)Cl",
    "toluene": "Cc1ccccc1",
    "water": "O",
    "o-xylene": "Cc1ccccc1C",
    "m-xylene": "Cc1cccc(C)c1",
    "p-xylene": "Cc1ccc(C)cc1",
    "xylene": "Cc1ccccc1C",  # defaults to o-xylene
    "2-heptanone": "CCCCCC(C)=O",
}


def _resolve_solvent(solvent: str) -> str:
    """
    Convert a solvent name or SMILES to SMILES.

    :param solvent: Solvent name (e.g., "ethanol") or SMILES (e.g., "CCO").
    :returns: SMILES string.
    :raises ValueError: If solvent name is not recognized and doesn't look like SMILES.
    """
    # Check if it's a known solvent name (case-insensitive)
    lower = solvent.lower().strip()
    if lower in COMMON_SOLVENTS:
        return COMMON_SOLVENTS[lower]

    # Check if it looks like a SMILES (contains uppercase letters typical of SMILES)
    # SMILES typically have uppercase C, N, O, S, F, Cl, Br, I, P, etc.
    if any(c in solvent for c in "CNOSFIBPcnos[]()=#"):
        return solvent

    # Unrecognized name
    raise ValueError(
        f"Unrecognized solvent '{solvent}'. "
        f"Use a common name (e.g., 'ethanol', 'water', 'thf') or provide a SMILES string. "
        f"Available names: {', '.join(sorted(COMMON_SOLVENTS.keys()))}"
    )


@dataclass(frozen=True, slots=True)
class SolubilityValue:
    """A solubility measurement at a specific temperature."""

    temperature: float
    """Temperature in Kelvin."""

    solubility: float
    """Solubility in log(mol/L)."""

    uncertainty: float | None = None
    """Uncertainty in the solubility prediction."""


@dataclass(frozen=True, slots=True)
class SolubilityEntry:
    """Solubility results for a single solvent."""

    solvent: str
    """Solvent SMILES."""

    values: tuple[SolubilityValue, ...]
    """Solubility values at each temperature."""


@register_result("solubility")
class SolubilityResult(WorkflowResult):
    """Result from an aqueous solubility workflow."""

    _stjames_class = stjames.SolubilityWorkflow

    def __repr__(self) -> str:
        solvents = [s.solvent for s in self.solubilities]
        return f"<SolubilityResult solvents={solvents}>"

    @property
    def solubilities(self) -> list[SolubilityEntry]:
        """Solubility results per solvent, with each value paired to its temperature."""
        temps = list(self._workflow.temperatures)
        entries = []
        for solvent, result in self._workflow.solubilities.items():
            values = tuple(
                SolubilityValue(
                    temperature=temps[i],
                    solubility=result.solubilities[i],
                    uncertainty=result.uncertainties[i] if result.uncertainties else None,
                )
                for i in range(len(temps))
            )
            entries.append(SolubilityEntry(solvent=solvent, values=values))
        return entries


def submit_solubility_workflow(
    initial_smiles: str | MoleculeInput,
    method: Literal["fastsolv", "kingfisher", "esol"] = "fastsolv",
    solvents: list[str] | None = None,
    temperatures: list[float] | None = None,
    name: str = "Solubility Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a solubility workflow to the API.

    :param initial_smiles: Molecule to calculate solubility for. Accepts a SMILES
        string or any molecule type (RowanMolecule, stjames.Molecule, RDKit Mol, or dict).
        The molecule must have a SMILES string associated with it, as solubility models
        are 2D/SMILES-based and do not use 3D coordinates.
    :param method: Solubility prediction method:
        - "fastsolv": ML-based solid solubility. Supports arbitrary solvents and temperatures.
        - "kingfisher": ML-based aqueous solubility. Water only, 298.15K only.
        - "esol": ESOL regression for aqueous solubility. Water only, 298.15K only.
    :param solvents: List of solvent names or SMILES. Common names like "ethanol",
        "water", "thf" are recognized (see COMMON_SOLVENTS). For fastsolv, any solvent
        SMILES is accepted. For kingfisher/esol, must be ["water"] or ["O"].
    :param temperatures: List of temperatures in Kelvin. For fastsolv, any temperatures.
        For kingfisher/esol, must be [298.15] (room temperature).
    :param name: Name of the workflow.
    :param folder_uuid: UUID of the folder to place the workflow in.
    :param max_credits: Maximum number of credits to use for the workflow.
    :returns: Workflow object representing the submitted workflow.
    :raises ValueError: If the molecule has no SMILES, or solvents/temperatures are
        incompatible with the method.
    :raises requests.HTTPError: If the request to the API fails.
    """
    initial_smiles = extract_smiles(initial_smiles)
    # Resolve solvent names to SMILES
    if solvents is not None:
        solvents = [_resolve_solvent(s) for s in solvents]

    # Method-specific defaults and validation
    match method:
        case "kingfisher" | "esol":
            if solvents is None:
                solvents = ["O"]
            elif solvents != ["O"]:
                raise ValueError(
                    f"Method '{method}' only supports aqueous solubility. "
                    f"solvents must be ['water'] or ['O'], got {solvents}"
                )
            if temperatures is None:
                temperatures = [298.15]
            elif len(temperatures) != 1 or abs(temperatures[0] - 298.15) > 0.1:
                raise ValueError(
                    f"Method '{method}' only supports room temperature (298.15K). "
                    f"Got {temperatures}"
                )
        case "fastsolv":
            if solvents is None:
                # Default: hexane, toluene, THF, ethyl acetate, ethanol, acetonitrile
                solvents = ["CCCCCC", "Cc1ccccc1", "C1CCCO1", "CC(=O)OCC", "CCO", "CC#N"]
            if temperatures is None:
                temperatures = [273.15, 298.15, 323.15, 348.15, 373.15]

    workflow = stjames.SolubilityWorkflow(
        initial_smiles=initial_smiles,
        solubility_method=method,
        solvents=solvents,
        temperatures=temperatures,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "solubility",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_smiles": initial_smiles,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())
