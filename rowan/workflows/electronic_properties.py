"""Electronic properties workflow - calculate electronic properties."""

from typing import Any

import stjames
from rdkit import Chem

from ..utils import api_client
from .base import RdkitMol, StJamesMolecule, Workflow, WorkflowResult, register_result


@register_result("electronic_properties")
class ElectronicPropertiesResult(WorkflowResult):
    """Result from an electronic properties workflow."""

    _stjames_class = stjames.ElectronicPropertiesWorkflow

    def __repr__(self) -> str:
        homo = getattr(self._workflow, "homo_energy", None)
        lumo = getattr(self._workflow, "lumo_energy", None)
        return f"<ElectronicPropertiesResult HOMO={homo} LUMO={lumo}>"


def submit_electronic_properties_workflow(
    initial_molecule: dict[str, Any] | StJamesMolecule | RdkitMol,
    method: stjames.Method | str = "b97_3c",
    basis_set: str | None = None,
    compute_density_cube: bool = True,
    compute_electrostatic_potential_cube: bool = True,
    compute_num_occupied_orbitals: int = 1,
    compute_num_virtual_orbitals: int = 1,
    name: str = "Electronic Properties Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits an electronic properties workflow to the API.

    :param initial_molecule: The molecule to calculate electronic properties for.
    :param method: The method to use for the calculation.
    :param basis_set: The basis set to use (if any).
    :param compute_density_cube: Whether to compute the density cube.
    :param compute_electrostatic_potential_cube: Whether to compute the electrostatic
        potential cube.
    :param compute_num_occupied_orbitals: Number of occupied orbitals to save.
    :param compute_num_virtual_orbitals: Number of virtual orbitals to save.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if isinstance(initial_molecule, StJamesMolecule):
        initial_molecule = initial_molecule.model_dump(mode="json")
    elif isinstance(initial_molecule, Chem.rdchem.Mol | Chem.rdchem.RWMol):
        initial_molecule = StJamesMolecule.from_rdkit(initial_molecule, cid=0).model_dump(
            mode="json"
        )

    if isinstance(method, str):
        method = stjames.Method(method)

    settings = stjames.Settings(method=method, basis_set=basis_set)

    workflow = stjames.ElectronicPropertiesWorkflow(
        initial_molecule=initial_molecule,
        settings=settings,
        compute_density_cube=compute_density_cube,
        compute_electrostatic_potential_cube=compute_electrostatic_potential_cube,
        compute_num_occupied_orbitals=compute_num_occupied_orbitals,
        compute_num_virtual_orbitals=compute_num_virtual_orbitals,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "electronic_properties",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": initial_molecule,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


__all__ = ["ElectronicPropertiesResult", "submit_electronic_properties_workflow"]
