"""Basic calculation workflow - perform quantum chemical calculations."""

from typing import Any

import stjames
from rdkit import Chem

from ..utils import api_client
from .base import RdkitMol, StJamesMolecule, Workflow, WorkflowResult, register_result


@register_result("basic_calculation")
class BasicCalculationResult(WorkflowResult):
    """Result from a basic calculation workflow."""

    # TODO: Add properties for energy, geometry, etc. once we figure out
    # how to fetch/embed calculation results from the calculation_uuid

    _stjames_class = stjames.BasicCalculationWorkflow

    def __repr__(self) -> str:
        calc_uuid = getattr(self._workflow, "calculation_uuid", None)
        return f"<BasicCalculationResult calculation_uuid={calc_uuid!r}>"


def submit_basic_calculation_workflow(
    initial_molecule: dict[str, Any] | StJamesMolecule | RdkitMol,
    method: stjames.Method | str = "uma_m_omol",
    basis_set: stjames.BasisSet | str | None = None,
    tasks: list[str] | None = None,
    mode: str = "auto",
    engine: str | None = None,
    name: str = "Basic Calculation Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submit a basic calculation workflow to the API.

    :param initial_molecule: The molecule to perform the calculation on.
    :param method: The method to use for the calculation.
    :param basis_set: The basis_set to use (if any).
    :param tasks: A list of tasks to perform for the calculation.
    :param mode: The mode to run the calculation in.
    :param engine: The engine to use for the calculation.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if not tasks:
        tasks = ["optimize"]

    if isinstance(initial_molecule, StJamesMolecule):
        initial_molecule = initial_molecule.model_dump(mode="json")
    elif isinstance(initial_molecule, Chem.rdchem.Mol | Chem.rdchem.RWMol):
        initial_molecule = StJamesMolecule.from_rdkit(initial_molecule, cid=0)

    if isinstance(method, str):
        method = stjames.Method(method)

    workflow = stjames.BasicCalculationWorkflow(
        initial_molecule=initial_molecule,
        settings=stjames.Settings(
            method=method,
            basis_set=basis_set,
            tasks=tasks,
            mode=mode,
        ),
        engine=engine or method.default_engine(),
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "basic_calculation",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": initial_molecule,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


__all__ = ["BasicCalculationResult", "submit_basic_calculation_workflow"]
