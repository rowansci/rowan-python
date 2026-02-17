"""pKa workflow - predict acid/base dissociation constants."""

from dataclasses import dataclass
from typing import Any, Literal

import stjames
from rdkit import Chem

from ..utils import api_client
from .base import RdkitMol, StJamesMolecule, Workflow, WorkflowResult, register_result


@dataclass(frozen=True, slots=True)
class pKaMicrostate:
    """A microstate from a pKa calculation."""

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
        return f"<pKaResult strongest_acid={acid} strongest_base={base}>"

    @property
    def strongest_acid(self) -> float | None:
        """Strongest acidic site pKa value."""
        return self._workflow.strongest_acid

    @property
    def strongest_base(self) -> float | None:
        """Strongest basic site pKa value."""
        return self._workflow.strongest_base

    @property
    def conjugate_acids(self) -> list[pKaMicrostate]:
        """List of conjugate acid microstates with pKa values."""
        return [
            pKaMicrostate(
                atom_index=m.atom_index,
                pka=m.pka,
                smiles=m.smiles,
                delta_g=m.deltaG,
                uncertainty=m.uncertainty,
            )
            for m in self._workflow.conjugate_acids
        ]

    @property
    def conjugate_bases(self) -> list[pKaMicrostate]:
        """List of conjugate base microstates with pKa values."""
        return [
            pKaMicrostate(
                atom_index=m.atom_index,
                pka=m.pka,
                smiles=m.smiles,
                delta_g=m.deltaG,
                uncertainty=m.uncertainty,
            )
            for m in self._workflow.conjugate_bases
        ]


def submit_pka_workflow(
    initial_molecule: dict[str, Any] | StJamesMolecule | RdkitMol | str,
    pka_range: tuple[int, int] = (2, 12),
    method: Literal["aimnet2_wagen2024", "chemprop_nevolianis2025"] = "aimnet2_wagen2024",
    solvent: str = "water",
    deprotonate_elements: list[int] | None = None,
    protonate_elements: list[int] | None = None,
    mode: str = "careful",
    name: str = "pKa Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a pKa workflow to the API.

    :param initial_molecule: The molecule to calculate the pKa of.
        Valid input types include `stjames.Molecule`, RDKit molecule, and SMILES.
    :param pka_range: The range of pKa values to calculate.
    :param method: The algorithm used to compute pKa values.
    :param solvent: The solvent in which pKa values will be computed.
    :param deprotonate_elements: The elements to deprotonate. Given by atomic number.
    :param protonate_elements: The elements to protonate. Given by atomic number.
    :param mode: The mode to run the calculation in.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    initial_smiles: str = ""
    initial_stjames_mol: StJamesMolecule | None = None

    if isinstance(initial_molecule, dict):
        initial_stjames_mol = StJamesMolecule.model_validate(initial_molecule)
    elif isinstance(initial_molecule, StJamesMolecule):
        initial_stjames_mol = initial_molecule
    elif isinstance(initial_molecule, Chem.rdchem.Mol | Chem.rdchem.RWMol):
        initial_stjames_mol = StJamesMolecule.from_rdkit(initial_molecule, cid=0)
    elif isinstance(initial_molecule, str):
        initial_smiles = initial_molecule

    protonate_elements = protonate_elements or [7]
    deprotonate_elements = deprotonate_elements or [7, 8, 16]

    workflow = stjames.pKaWorkflow(
        initial_molecule=initial_stjames_mol,
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
        "initial_molecule": initial_stjames_mol.model_dump(mode="json")
        if initial_stjames_mol
        else None,
        "initial_smiles": initial_smiles,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


__all__ = ["pKaMicrostate", "pKaResult", "submit_pka_workflow"]
