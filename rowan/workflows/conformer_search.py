"""Conformer search workflow - find low-energy molecular conformations."""

from typing import Any

import stjames
from rdkit import Chem

from ..utils import api_client
from .base import RdkitMol, StJamesMolecule, Workflow, WorkflowResult, register_result


@register_result("conformer_search")
class ConformerSearchResult(WorkflowResult):
    """Result from a conformer search workflow."""

    _stjames_class = stjames.ConformerSearchWorkflow

    def __repr__(self) -> str:
        energies = self.energies
        n = len(energies)
        e_min = min(energies) if energies else None
        e_max = max(energies) if energies else None
        return f"<ConformerSearchResult conformers={n} energy_range=({e_min}, {e_max})>"

    @property
    def conformer_uuids(self) -> list[list[str | None]]:
        """List of conformer UUIDs (nested for multistage optimization)."""
        return self._workflow.conformer_uuids

    @property
    def energies(self) -> list[float]:
        """Conformer energies."""
        return list(self._workflow.energies)


def submit_conformer_search_workflow(
    initial_molecule: dict[str, Any] | StJamesMolecule | RdkitMol,
    conf_gen_settings: stjames.ConformerGenSettings,
    final_method: stjames.Method | str = "aimnet2_wb97md3",
    solvent: str | None = None,
    transition_state: bool = False,
    name: str = "Conformer Search Workflow",
    folder_uuid: str | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a conformer search workflow to the API.

    :param initial_molecule: The molecule to perform the conformer search on.
    :param conf_gen_settings: settings for conformer generation
    :param final_method: The method to use for the final optimization.
    :param solvent: The solvent to use for the final optimization.
    :param transition_state: Whether to optimize the transition state.
    :param name: The name of the workflow.
    :param folder_uuid: The UUID of the folder to place the workflow in.
    :param max_credits: The maximum number of credits to use for the workflow.
    :return: A Workflow object representing the submitted workflow.
    :raises requests.HTTPError: if the request to the API fails.
    """
    if isinstance(initial_molecule, StJamesMolecule):
        initial_molecule = initial_molecule.model_dump(mode="json")
    elif isinstance(initial_molecule, Chem.rdchem.Mol | Chem.rdchem.RWMol):
        initial_molecule = StJamesMolecule.from_rdkit(initial_molecule, cid=0)

    if isinstance(final_method, str):
        final_method = stjames.Method(final_method)

    solvent_model = None
    if solvent:
        solvent_model = "alpb" if final_method in stjames.XTB_METHODS else "cpcm"

    opt_settings = stjames.Settings(
        method=final_method,
        tasks=["optimize"],
        mode=stjames.Mode.AUTO,
        solvent_settings={"solvent": solvent, "model": solvent_model} if solvent else None,
        opt_settings={"transition_state": transition_state, "constraints": []},
    )

    msos = stjames.MultiStageOptSettings(
        mode=stjames.Mode.MANUAL,
        xtb_preopt=True,
        optimization_settings=[opt_settings],
    )

    workflow = stjames.ConformerSearchWorkflow(
        initial_molecule=initial_molecule,
        multistage_opt_settings=msos,
        conf_gen_settings=conf_gen_settings,
        solvent=solvent,
        transition_state=transition_state,
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "conformer_search",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": initial_molecule,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())


__all__ = ["ConformerSearchResult", "submit_conformer_search_workflow"]
