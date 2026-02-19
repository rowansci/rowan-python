"""Conformer search workflow - find low-energy molecular conformations."""

from typing import Any

import stjames
from rdkit import Chem

from ..calculation import Calculation, retrieve_calculation
from ..molecule import Molecule
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
        """Conformer energies (Hartree)."""
        return list(self._workflow.energies)

    @property
    def radii_of_gyration(self) -> list[float]:
        """Radius of gyration for each conformer (Å)."""
        return [p.radius_of_gyration for p in self._workflow.conformer_properties]

    @property
    def sasa(self) -> list[float]:
        """Solvent accessible surface area for each conformer (Ų)."""
        return [p.solvent_accessible_surface_area for p in self._workflow.conformer_properties]

    @property
    def polar_sasa(self) -> list[float]:
        """Polar solvent accessible surface area for each conformer (Ų)."""
        return [
            p.polar_solvent_accessible_surface_area for p in self._workflow.conformer_properties
        ]

    def get_conformers(self, n: int | None = None) -> list[Molecule]:
        """
        Fetch conformer molecules.

        :param n: Number of conformers to fetch (default: all). Conformers are
            ordered by energy, so n=5 returns the 5 lowest-energy conformers.
        :return: List of Molecule objects.

        Note: Makes one API call per conformer.
        """
        count = len(self.conformer_uuids) if n is None else min(n, len(self.conformer_uuids))
        molecules = []
        for i in range(count):
            calc = self.get_conformer(i)
            if calc.molecule:
                molecules.append(calc.molecule)
        return molecules

    def get_conformer(self, index: int, stage: int = -1) -> Calculation:
        """
        Fetch a conformer's calculation data by index.

        :param index: The conformer index (0-based).
        :param stage: The optimization stage (-1 for final stage).
        :return: A Calculation object with molecule and energy data.
        :raises IndexError: If the index is out of range.
        :raises ValueError: If the conformer UUID is None.
        """
        uuids = self.conformer_uuids
        if index < 0 or index >= len(uuids):
            raise IndexError(f"Conformer index {index} out of range (0-{len(uuids) - 1})")

        stage_uuids = uuids[index]
        uuid = stage_uuids[stage]
        if uuid is None:
            raise ValueError(f"Conformer {index} has no calculation at stage {stage}")

        return retrieve_calculation(uuid)


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
