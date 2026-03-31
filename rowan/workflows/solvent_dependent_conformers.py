"""Solvent-dependent conformers workflow - conformer search with multi-solvent scoring."""

from dataclasses import dataclass

import stjames
from stjames import ConformerGenSettingsUnion, iMTDSettings

from ..folder import Folder
from ..types import MoleculeInput
from ..utils import api_client
from .base import Solvent, Workflow, WorkflowResult, molecule_to_dict, register_result

_DEFAULT_SOLVENTS: list[Solvent] = [
    Solvent.HEXANE,
    Solvent.OCTANOL,
    Solvent.CHLOROFORM,
    Solvent.DIMETHYLSULFOXIDE,
    Solvent.WATER,
]


@dataclass(frozen=True, slots=True)
class SolventDependentConformerProperties:
    """Conformer ensemble properties for a single solvent.

    :param solvent_accessible_surface_area: Average SASA (A^2).
    :param polar_solvent_accessible_surface_area: Average polar SASA for non-C/H atoms (A^2).
    :param radius_of_gyration: Radius of gyration (A).
    """

    solvent_accessible_surface_area: float
    polar_solvent_accessible_surface_area: float
    radius_of_gyration: float


@dataclass(frozen=True, slots=True)
class SolventDependentConformer:
    """A single conformer scored across multiple solvents.

    :param calculation_uuid: UUID of the underlying calculation.
    :param free_energy_by_solvent: Absolute free energy per solvent (Hartree).
    :param relative_free_energy_by_solvent: Free energy relative to lowest conformer per solvent
        (kcal/mol).
    :param population_by_solvent: Boltzmann population per solvent (0-1).
    """

    calculation_uuid: str
    free_energy_by_solvent: dict[Solvent, float]
    relative_free_energy_by_solvent: dict[Solvent, float]
    population_by_solvent: dict[Solvent, float]


@register_result("solvent_dependent_conformers")
class SolventDependentConformersResult(WorkflowResult):
    """Result from a solvent-dependent conformers workflow."""

    _stjames_class = stjames.SolventDependentConformersWorkflow

    def __repr__(self) -> str:
        return (
            f"<SolventDependentConformersResult conformers={self.num_conformers}"
            f" solvents={len(self.solvents)}>"
        )

    @property
    def num_conformers(self) -> int:
        """Number of conformers found."""
        return len(self._workflow.conformers)

    @property
    def solvents(self) -> list[Solvent]:
        """Solvents used for scoring."""
        return list(self._workflow.solvents)

    @property
    def conformers(self) -> list[SolventDependentConformer]:
        """Conformers with per-solvent energies and populations."""
        return [
            SolventDependentConformer(
                calculation_uuid=str(c.calculation),
                free_energy_by_solvent=dict(c.free_energy_by_solvent),
                relative_free_energy_by_solvent=dict(c.relative_free_energy_by_solvent),
                population_by_solvent=dict(c.population_by_solvent),
            )
            for c in self._workflow.conformers
        ]

    @property
    def per_solvent_properties(self) -> dict[Solvent, SolventDependentConformerProperties]:
        """Aggregate ensemble properties per solvent (SASA, polar SASA, radius of gyration)."""
        return {
            solvent: SolventDependentConformerProperties(
                solvent_accessible_surface_area=p.solvent_accessible_surface_area,
                polar_solvent_accessible_surface_area=p.polar_solvent_accessible_surface_area,
                radius_of_gyration=p.radius_of_gyration,
            )
            for solvent, p in self._workflow.per_solvent_properties.items()
        }

    @property
    def relative_free_energy_by_solvent(self) -> dict[Solvent, float]:
        """Relative transfer free energy by solvent (kcal/mol)."""
        return dict(self._workflow.relative_free_energy_by_solvent)


def submit_solvent_dependent_conformers_workflow(
    initial_molecule: MoleculeInput,
    solvents: list[Solvent] | None = None,
    conf_gen_settings: ConformerGenSettingsUnion | None = None,
    energy_window: float = 30,
    name: str = "Solvent-Dependent Conformers",
    folder_uuid: str | None = None,
    folder: Folder | None = None,
    max_credits: int | None = None,
    webhook_url: str | None = None,
) -> Workflow:
    """
    Submits a solvent-dependent conformers workflow to the API.

    Generates conformers and scores them across multiple solvents using CPCM-X,
    enabling prediction of solvent-dependent conformer populations and transfer
    free energies. The input molecule must have 3D coordinates.

    :param initial_molecule: Molecule with 3D coordinates to perform the conformer search on.
    :param solvents: Solvents to score conformers in (``rowan.Solvent`` enum).
        Defaults to hexane, octanol, chloroform, DMSO, and water.
    :param conf_gen_settings: Conformer generation settings. Defaults to
        ``rowan.iMTDSettings`` with GFN-FF and 30 kcal/mol energy window.
        Other options: ``rowan.ETKDGSettings``, ``rowan.iMTDGCSettings``,
        ``rowan.LyrebirdSettings``.
    :param energy_window: Energy window for conformer generation (kcal/mol).
        Only used when ``conf_gen_settings`` is None.
    :param name: Name of the workflow.
    :param folder_uuid: UUID of the folder to place the workflow in.
    :param folder: Folder object to store the workflow in.
    :param max_credits: Maximum number of credits to use for the workflow.
    :param webhook_url: URL that Rowan will POST to when the workflow completes.
    :returns: Workflow object representing the submitted workflow.
    :raises ValueError: If both folder and folder_uuid are provided.
    :raises requests.HTTPError: If the request to the API fails.
    """
    if folder and folder_uuid:
        raise ValueError("Provide either `folder` or `folder_uuid`, not both.")
    if folder:
        folder_uuid = folder.uuid

    if conf_gen_settings is None:
        conf_gen_settings = iMTDSettings(
            max_confs=None,
            speed="normal",
            solvent_settings={"solvent": "water", "model": "alpb"},
            reopt=False,
            mtd_method="gfn_ff",
            energy_window=energy_window,
        )

    mol_dict = molecule_to_dict(initial_molecule)

    workflow = stjames.SolventDependentConformersWorkflow(
        initial_molecule=mol_dict,
        conf_gen_settings=conf_gen_settings,
        solvents=solvents if solvents is not None else _DEFAULT_SOLVENTS,
    )

    data = {
        "workflow_type": "solvent_dependent_conformers",
        "workflow_data": workflow.model_dump(mode="json"),
        "initial_molecule": mol_dict,
        "name": name,
        "folder_uuid": folder_uuid,
        "max_credits": max_credits,
        "webhook_url": webhook_url,
    }

    with api_client() as client:
        response = client.post("/workflow", json=data)
        response.raise_for_status()
        return Workflow(**response.json())
