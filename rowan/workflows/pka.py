"""pKa workflow - predict acid/base dissociation constants."""

from dataclasses import dataclass
from typing import Literal

import stjames

from ..calculation import Calculation, retrieve_calculation
from ..utils import api_client
from .base import (
    Mode,
    MoleculeInput,
    SolventInput,
    Workflow,
    WorkflowResult,
    molecule_to_dict,
    register_result,
)


@dataclass(frozen=True, slots=True)
class pKaMicrostate:
    """
    A microstate from a pKa calculation.

    Note: Available fields depend on the pKa method used.
    - aimnet2_wagen2024: has delta_g
    - chemprop_nevolianis2025: has smiles and uncertainty
    """

    atom_index: int
    pka: float
    smiles: str | None = None
    delta_g: float | None = None
    uncertainty: float | None = None


@register_result("pka")
class pKaResult(WorkflowResult):
    """Result from a pKa workflow."""

    _stjames_class = stjames.pKaWorkflow

    def __repr__(self) -> str:
        acid = self.strongest_acid
        base = self.strongest_base
        parts = []
        if acid is not None:
            parts.append(f"acid={acid:.2f}")
        if base is not None:
            parts.append(f"base={base:.2f}")
        return f"<pKaResult {' '.join(parts)}>" if parts else "<pKaResult>"

    @property
    def strongest_acid(self) -> float | None:
        """Strongest acidic site pKa value."""
        return getattr(self._workflow, "strongest_acid", None)

    @property
    def strongest_base(self) -> float | None:
        """Strongest basic site pKa value."""
        return getattr(self._workflow, "strongest_base", None)

    @property
    def conjugate_acids(self) -> list[pKaMicrostate]:
        """List of conjugate acid microstates with pKa values."""
        return [self._make_microstate(m) for m in self._workflow.conjugate_acids]

    @property
    def conjugate_bases(self) -> list[pKaMicrostate]:
        """List of conjugate base microstates with pKa values."""
        return [self._make_microstate(m) for m in self._workflow.conjugate_bases]

    def _make_microstate(self, m) -> pKaMicrostate:
        """Convert a stjames microstate to pKaMicrostate."""
        return pKaMicrostate(
            atom_index=m.atom_index,
            pka=m.pka,
            smiles=getattr(m, "smiles", None),
            delta_g=getattr(m, "deltaG", None),
            uncertainty=getattr(m, "uncertainty", None),
        )

    def _get_structure_uuids(self) -> list[str]:
        """Extract structure UUIDs from workflow data."""
        structures = getattr(self._workflow, "structures", None) or []
        uuids = []
        for s in structures:
            if s:
                uuid = s.get("uuid") if isinstance(s, dict) else getattr(s, "uuid", None)
                if uuid:
                    uuids.append(uuid)
        return uuids

    @property
    def structures(self) -> list[Calculation]:
        """
        Optimized structure calculations (lazily fetched).

        Only available for aimnet2_wagen2024 method (3D structure-based).

        Note: Makes one API call per structure on first access.
        Results are cached. Call clear_cache() to refresh.

        :raises ValueError: If method is chemprop_nevolianis2025 (no structures).
        """
        method = getattr(self._workflow, "microscopic_pka_method", None)
        if method == "chemprop_nevolianis2025":
            raise ValueError(
                "chemprop_nevolianis2025 is SMILES-based and does not produce structures. "
                "Use aimnet2_wagen2024 for 3D structure calculations."
            )
        if "structures" not in self._cache:
            uuids = self._get_structure_uuids()
            self._cache["structures"] = [retrieve_calculation(uuid) for uuid in uuids]
        return self._cache["structures"]


def submit_pka_workflow(
    initial_molecule: MoleculeInput | str,
    pka_range: tuple[int, int] = (2, 12),
    method: Literal["aimnet2_wagen2024", "chemprop_nevolianis2025"] = "aimnet2_wagen2024",
    solvent: SolventInput = "water",
    deprotonate_elements: list[int] | None = None,
    protonate_elements: list[int] | None = None,
    mode: Mode = Mode.CAREFUL,
    name: str = "pKa Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a pKa workflow to the API.

    :param initial_molecule: The molecule to calculate the pKa of.
        Accepts Molecule, stjames.Molecule, RDKit Mol, dict, or SMILES string.
    :param pka_range: The range of pKa values to calculate.
    :param method: The algorithm used to compute pKa values.
    :param solvent: The solvent in which pKa values will be computed.
    :param deprotonate_elements: Elements to deprotonate (atomic numbers).
    :param protonate_elements: Elements to protonate (atomic numbers).
    :param mode: The mode to run the calculation in.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises ValueError: If method and input type don't match.
    :raises requests.HTTPError: if the request to the API fails.
    """
    initial_smiles: str = ""
    mol_dict: dict | None = None

    if isinstance(initial_molecule, str):
        initial_smiles = initial_molecule
    else:
        mol_dict = molecule_to_dict(initial_molecule)

    # Validate method/input compatibility
    if method == "aimnet2_wagen2024" and not mol_dict:
        raise ValueError(
            "aimnet2_wagen2024 requires a 3D structure. Provide a Molecule, not a SMILES string."
        )
    if method == "chemprop_nevolianis2025" and not initial_smiles:
        raise ValueError(
            "chemprop_nevolianis2025 requires a SMILES string. "
            "Provide a SMILES string, not a 3D structure."
        )

    protonate_elements = protonate_elements or [7]
    deprotonate_elements = deprotonate_elements or [7, 8, 16]

    workflow = stjames.pKaWorkflow(
        initial_molecule=mol_dict,
        initial_smiles=initial_smiles,
        pka_range=pka_range,
        deprotonate_elements=deprotonate_elements,
        protonate_elements=protonate_elements,
        mode=mode,
        solvent=solvent,
        microscopic_pka_method=method,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "pka",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": mol_dict,
        "initial_smiles": initial_smiles,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


__all__ = ["pKaMicrostate", "pKaResult", "submit_pka_workflow"]
