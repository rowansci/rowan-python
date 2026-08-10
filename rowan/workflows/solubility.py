"""Solubility workflow - predict molecular solubility in various solvents."""

from dataclasses import dataclass
from typing import Literal

import stjames
from rdkit import Chem

from ..folder import Folder
from ..utils import api_client
from .base import (
    SMILES,
    Workflow,
    WorkflowResult,
    extract_smiles,
    register_result,
    submit_workflow_group,
)

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
    :raises ValueError: If solvent is not a recognized name or valid SMILES.
    """
    # Check if it's a known solvent name (case-insensitive)
    lower = solvent.lower().strip()
    if lower in COMMON_SOLVENTS:
        return COMMON_SOLVENTS[lower]

    # Validate as SMILES using RDKit
    if Chem.MolFromSmiles(solvent) is not None:
        return solvent

    raise ValueError(
        f"Unrecognized solvent '{solvent}'. "
        f"Use a common name (e.g., 'ethanol', 'water', 'thf') or provide a valid SMILES string. "
        f"Available names: {', '.join(sorted(COMMON_SOLVENTS.keys()))}"
    )


@dataclass(frozen=True, slots=True)
class SolubilityValue:
    """
    Solubility measurement at a specific temperature.

    :param temperature: Temperature in Kelvin.
    :param solubility: Solubility in log(mol/L).
    :param uncertainty: Uncertainty in the solubility prediction.
    """

    temperature: float
    solubility: float
    uncertainty: float | None = None


@dataclass(frozen=True, slots=True)
class SolubilityEntry:
    """
    Solubility results for a single solvent.

    :param solvent: Solvent SMILES.
    :param values: Solubility values at each temperature.
    """

    solvent: str
    values: tuple[SolubilityValue, ...]


@register_result("solubility")
class SolubilityResult(WorkflowResult):
    """Result from an aqueous-solubility workflow."""

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


def _prepare_solubility_settings(
    method: Literal["fastsolv", "kingfisher", "esol"],
    solvents: list[str] | None,
    temperatures: list[float] | None,
) -> tuple[list[str], list[float]]:
    """Resolve and validate method-specific solubility settings."""
    resolved_solvents = (
        [_resolve_solvent(solvent) for solvent in solvents] if solvents is not None else None
    )

    match method:
        case "kingfisher" | "esol":
            if resolved_solvents is None:
                resolved_solvents = ["O"]
            elif resolved_solvents != ["O"]:
                raise ValueError(
                    f"Method '{method}' only supports aqueous solubility. "
                    f"solvents must be ['water'] or ['O'], got {resolved_solvents}"
                )
            if temperatures is None:
                temperatures = [298.15]
            elif len(temperatures) != 1 or abs(temperatures[0] - 298.15) > 0.1:
                raise ValueError(
                    f"Method '{method}' only supports room temperature (298.15K). "
                    f"Got {temperatures}"
                )
        case "fastsolv":
            if resolved_solvents is None:
                resolved_solvents = [
                    "CCCCCC",
                    "Cc1ccccc1",
                    "C1CCCO1",
                    "CC(=O)OCC",
                    "CCO",
                    "CC#N",
                ]
            if temperatures is None:
                temperatures = [273.15, 298.15, 323.15, 348.15, 373.15]

    return resolved_solvents, temperatures


def submit_solubility_workflow(
    initial_smiles: SMILES,
    method: Literal["fastsolv", "kingfisher", "esol"] = "fastsolv",
    solvents: list[str] | None = None,
    temperatures: list[float] | None = None,
    name: str = "Solubility Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
    is_draft: bool = False,
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
    :param folder: Folder object to store the workflow in.
    :param max_credits: Maximum number of credits to use for the workflow.
    :param webhook_url: URL that Rowan will POST to when the workflow completes.
    :param is_draft: If True, submit the workflow as a draft without starting execution.
    :returns: Workflow object representing the submitted workflow.
    :raises ValueError: If the molecule has no SMILES, or solvents/temperatures are
        incompatible with the method.
    :raises requests.HTTPError: If the request to the API fails.
    """
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder:
        folder_uuid = folder.uuid
    initial_smiles = extract_smiles(initial_smiles)
    solvents, temperatures = _prepare_solubility_settings(method, solvents, temperatures)

    workflow = stjames.SolubilityWorkflow(
        initial_smiles=initial_smiles,
        solubility_method=method,
        solvents=solvents,
        temperatures=temperatures,
    )

    data = {
        "workflow_type": "solubility",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_smiles": initial_smiles,
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


def submit_solubility_workflow_group(
    initial_smileses: list[SMILES],
    method: Literal["fastsolv", "kingfisher", "esol"] = "fastsolv",
    solvents: list[str] | None = None,
    temperatures: list[float] | None = None,
    names: list[str] | None = None,
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
) -> list[Workflow]:
    """Submit a batch of solubility workflows as one submission group.

    All molecules use the same method, solvents, and temperatures. A batch may contain
    up to 5,000 molecules. Each molecule is represented by its own ``Workflow``, and all
    returned workflows share a ``submission_group_uuid``. Use ``batch_poll_status()`` to
    monitor their UUIDs together and ``retrieve_workflows()`` to retrieve their records.

    :param initial_smileses: nonempty list of up to 5,000 solute SMILES strings
    :param method: solubility prediction method
    :param solvents: solvent names or SMILES strings
    :param temperatures: temperatures in Kelvin
    :param names: optional workflow names; when provided, one per SMILES string
    :param folder_uuid: UUID of the folder in which to store the workflows
    :param folder: folder in which to store the workflows
    :param max_credits: maximum credits to use per workflow
    :param webhook_url: URL Rowan will POST to when each workflow completes
    :returns: submitted workflows in one submission group
    """
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder:
        folder_uuid = folder.uuid
    if not initial_smileses:
        raise ValueError("Provide at least one initial SMILES string.")

    solvents, temperatures = _prepare_solubility_settings(method, solvents, temperatures)
    workflow = stjames.SolubilityWorkflow(
        initial_smiles=initial_smileses[0],
        solubility_method=method,
        solvents=solvents,
        temperatures=temperatures,
    )
    workflow_data = workflow.model_dump(
        mode="json",
        exclude={"initial_smiles", "messages", "solubilities"},
    )

    return submit_workflow_group(
        workflow_type="solubility",
        workflow_data=workflow_data,
        initial_smileses=initial_smileses,
        names=names,
        folder_uuid=folder_uuid,
        max_credits=max_credits,
        webhook_url=webhook_url,
    )
