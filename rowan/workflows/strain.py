"""Strain workflow - calculate molecular strain energy."""

import math

import stjames

from ..calculation import Calculation, retrieve_calculation
from ..folder import Folder
from ..molecule import Molecule
from ..utils import api_client
from .base import (
    Message,
    MoleculeInput,
    Workflow,
    WorkflowResult,
    molecule_to_dict,
    parse_messages,
    register_result,
)
from .constants import BOLTZMANN_HARTREE_PER_K


@register_result("strain")
class StrainResult(WorkflowResult):
    """Result from a strain workflow."""

    _stjames_class = stjames.StrainWorkflow

    def __repr__(self) -> str:
        return f"<StrainResult strain={self.strain} kcal/mol>"

    @property
    def strain(self) -> float | None:
        """Computed strain energy (kcal/mol)."""
        return self._workflow.strain

    @property
    def conformer_uuids(self) -> list[str | None]:
        """UUIDs of conformer calculations."""
        return self._workflow.conformers

    @property
    def constrained_optimization_uuid(self) -> str | None:
        """UUID of the constrained optimization calculation."""
        return getattr(self._workflow, "constrained_optimization", None)

    @property
    def messages(self) -> list[Message]:
        """Any messages or warnings from the workflow."""
        return parse_messages(getattr(self._workflow, "messages", None))

    @property
    def constrained_optimization(self) -> Calculation | None:
        """The constrained optimization calculation."""
        if not (uuid := self.constrained_optimization_uuid):
            return None

        if "constrained_opt" not in self._cache:
            self._cache["constrained_opt"] = retrieve_calculation(uuid)
        return self._cache["constrained_opt"]

    @property
    def conformers(self) -> list[Calculation]:
        """All conformer calculations.

        .. note::
            Makes one API call per conformer on first access.
            Results are cached. Call clear_cache() to refresh.
        """
        if "all_conformers" not in self._cache:
            calcs = []
            for uuid in self.conformer_uuids:
                if uuid:
                    calcs.append(retrieve_calculation(uuid))
            self._cache["all_conformers"] = calcs
        return self._cache["all_conformers"]

    @property
    def conformer_energies(self) -> list[float | None]:
        """Energies for all conformers (Hartree)."""
        return [c.energy for c in self.conformers]

    @property
    def conformer_molecules(self) -> list[Molecule]:
        """Molecule objects for all conformers."""
        return [c.molecule for c in self.conformers if c.molecule]

    def get_boltzmann_weights(self, temperature: float = 300.0) -> list[float]:
        """
        Compute Boltzmann weights for conformers.

        :param temperature: Temperature in Kelvin (default: 300K).
        :returns: List of weights (sum to 1.0), excluding failed conformers.
        """
        energies = self.conformer_energies
        valid_energies = [e for e in energies if e is not None]

        if not valid_energies:
            return []

        min_e = min(valid_energies)
        weights = []
        for e in valid_energies:
            rel_e = e - min_e
            weights.append(math.exp(-rel_e / (BOLTZMANN_HARTREE_PER_K * temperature)))

        total = sum(weights)
        return [w / total for w in weights]


def submit_strain_workflow(
    initial_molecule: MoleculeInput,
    harmonic_constraint_spring_constant: float = 5.0,
    constrain_hydrogens: bool = False,
    conf_gen_settings: stjames.ConformerGenSettingsUnion | None = None,
    name: str = "Strain Workflow",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
) -> Workflow:
    """
    Submits a strain workflow to the API.

    :param initial_molecule: Molecule to calculate strain for.
    :param harmonic_constraint_spring_constant: Spring constant for harmonic
        constraints (kcal/mol/A). Default 5.0.
    :param constrain_hydrogens: Whether to constrain hydrogen positions. Default False.
    :param conf_gen_settings: Conformer generation settings. Defaults to ETKDG with
        max 50 conformers.
    :param name: Name of the workflow.
    :param folder_uuid: UUID of the folder to store the workflow in.
    :param folder: Folder object to store the workflow in.
    :param max_credits: Maximum number of credits to use for the workflow.
    :returns: Workflow object representing the submitted workflow.
    :raises requests.HTTPError: If the request to the API fails.
    """
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder:
        folder_uuid = folder.uuid
    initial_molecule = molecule_to_dict(initial_molecule)

    workflow = stjames.StrainWorkflow(
        initial_molecule=initial_molecule,
        harmonic_constraint_spring_constant=harmonic_constraint_spring_constant,
        constrain_hydrogens=constrain_hydrogens,
        conf_gen_settings=conf_gen_settings or stjames.ETKDGSettings(max_confs=50),
    )

    data = {
        "name": name,
        "folder_uuid": folder_uuid,
        "workflow_type": "strain",
        "workflow_data": workflow.model_dump(serialize_as_any=True, mode="json"),
        "initial_molecule": initial_molecule,
        "max_credits": max_credits,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())
