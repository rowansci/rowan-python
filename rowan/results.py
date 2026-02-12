"""
Result types for Rowan workflows.

This module provides typed result classes that wrap workflow output data,
giving users autocomplete and type safety when accessing results.

Design principles:
1. Flatten outputs - property-wise access, not nested object diving
2. Python native types - return float, str, list, not stjames objects
3. Curated properties - expose only what users need, hide internal complexity
4. Escape hatches - .data for raw dict access when needed
5. stjames under the hood - use for validation/parsing, don't expose to users
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, ClassVar

import stjames


@dataclass(frozen=True, slots=True)
class BDEEntry:
    """A bond dissociation energy result."""

    fragment_idxs: tuple[int, ...]
    energy: float | None = None


@dataclass(frozen=True, slots=True)
class CofoldingResult:
    """A protein cofolding result."""

    plddt: float | None = None
    ptm: float | None = None
    iptm: float | None = None


@dataclass(frozen=True, slots=True)
class DockingScore:
    """A docking pose with its score."""

    score: float
    pose: str | None = None
    complex_pdb: str | None = None
    posebusters_valid: bool = False
    strain: float | None = None
    rmsd: float | None = None


@dataclass(frozen=True, slots=True)
class MacropKaMicrostate:
    """A microstate from a macropKa calculation."""

    smiles: str
    energy: float
    charge: int


@dataclass(frozen=True, slots=True)
class MacropKaValue:
    """A macroscopic pKa value."""

    initial_charge: int
    final_charge: int
    pka: float


@dataclass(frozen=True, slots=True)
class NMRPeak:
    """An NMR peak."""

    nucleus: int
    shift: float
    atom_indices: tuple[int, ...]


@dataclass(frozen=True, slots=True)
class pKaMicrostate:
    """A microstate from a pKa calculation."""

    atom_index: int
    pka: float
    smiles: str | None = None
    delta_g: float | None = None
    uncertainty: float | None = None


@dataclass(frozen=True, slots=True)
class ProteinBinder:
    """A generated protein binder design."""

    sequence: str
    plddt: float | None = None
    ptm: float | None = None
    iptm: float | None = None


@dataclass(frozen=True, slots=True)
class SolubilityEntry:
    """Solubility results for a single solvent."""

    solvent: str
    solubilities: tuple[float, ...]
    uncertainties: tuple[float | None, ...]


@dataclass(frozen=True, slots=True)
class SpinState:
    """A spin state result."""

    multiplicity: int
    energy: float


@dataclass(frozen=True, slots=True)
class Tautomer:
    """A tautomer result."""

    energy: float
    weight: float | None = None
    predicted_relative_energy: float | None = None


@dataclass(slots=True)
class WorkflowResult:
    """
    Base class for workflow results.

    Wraps the raw workflow data dict and parses it into a stjames object
    for typed access to nested data.

    :param workflow_data: The raw data dict from the workflow
    :param workflow_type: The workflow type string
    :param workflow_metadata: Metadata dict (elapsed, credits_charged, timestamps)
    """

    workflow_data: dict
    workflow_type: str
    workflow_metadata: dict = field(default_factory=dict)
    _workflow: Any = field(default=None, init=False, repr=False)

    _stjames_class: ClassVar[type | None] = None

    def __post_init__(self) -> None:
        """Parse workflow data into stjames object for typed access."""
        if self._stjames_class is not None:
            obj = self._stjames_class.model_validate(self.workflow_data)  # type: ignore[attr-defined]
            object.__setattr__(self, "_workflow", obj)

    @property
    def data(self) -> dict:
        """Raw workflow data dict for fallback access."""
        return self.workflow_data

    @property
    def metadata(self) -> dict:
        """
        Workflow metadata: elapsed, credits_charged, timestamps.

        Example keys:
        - elapsed: runtime in seconds
        - credits_charged: credits used
        - created_at: workflow creation time
        - started_at: computation start time
        - completed_at: computation end time
        """
        return self.workflow_metadata


class ADMETResult(WorkflowResult):
    """Result from an ADMET workflow."""

    _stjames_class = stjames.ADMETWorkflow

    @property
    def properties(self) -> dict[str, float | int] | None:
        """ADMET properties (molecular weight, logP, TPSA, etc.)."""
        return self._workflow.properties


class AnalogueDockingResult(WorkflowResult):
    """Result from an analogue docking workflow."""

    _stjames_class = stjames.AnalogueDockingWorkflow

    @property
    def analogue_scores(self) -> dict[str, list[DockingScore]]:
        """Docking scores for each analogue SMILES."""
        return {
            smiles: [
                DockingScore(
                    score=s.score,
                    pose=s.pose,
                    complex_pdb=s.complex_pdb,
                    posebusters_valid=s.posebusters_valid,
                    strain=s.strain,
                    rmsd=s.rmsd,
                )
                for s in scores
            ]
            for smiles, scores in self._workflow.analogue_scores.items()
        }


class BasicCalculationResult(WorkflowResult):
    """Result from a basic calculation workflow."""

    # TODO: Add properties for energy, geometry, etc. once we figure out
    # how to fetch/embed calculation results from the calculation_uuid

    _stjames_class = stjames.BasicCalculationWorkflow


class BatchDockingResult(WorkflowResult):
    """Result from a batch docking workflow."""

    _stjames_class = stjames.BatchDockingWorkflow

    @property
    def best_scores(self) -> list[float | None]:
        """Best docking scores for each ligand."""
        return list(self._workflow.best_scores)


class BDEResult(WorkflowResult):
    """Result from a Bond-Dissociation Energy (BDE) workflow."""

    _stjames_class = stjames.BDEWorkflow

    @property
    def optimization_energy(self) -> float | None:
        """Energy of optimized initial molecule."""
        return self._workflow.optimization_energy

    @property
    def bdes(self) -> list[BDEEntry]:
        """Bond dissociation energies."""
        return [
            BDEEntry(
                fragment_idxs=tuple(b.fragment_idxs),
                energy=b.energy,
            )
            for b in self._workflow.bdes
        ]

    @property
    def energies(self) -> tuple[float | None, ...]:
        """Convenience: extract energies from bdes."""
        return self._workflow.energies


class ConformerSearchResult(WorkflowResult):
    """Result from a conformer search workflow."""

    _stjames_class = stjames.ConformerSearchWorkflow

    @property
    def conformer_uuids(self) -> list[list[str | None]]:
        """List of conformer UUIDs (nested for multistage optimization)."""
        return self._workflow.conformer_uuids

    @property
    def energies(self) -> list[float]:
        """Conformer energies."""
        return list(self._workflow.energies)


class DescriptorsResult(WorkflowResult):
    """Result from a descriptors workflow."""

    _stjames_class = stjames.DescriptorsWorkflow

    @property
    def descriptors(self) -> dict | None:
        """Computed molecular descriptors."""
        return self._workflow.descriptors


class DockingResult(WorkflowResult):
    """Result from a docking workflow."""

    _stjames_class = stjames.DockingWorkflow

    @property
    def scores(self) -> list[DockingScore]:
        """List of docking scores with poses."""
        return [
            DockingScore(
                score=s.score,
                pose=s.pose,
                complex_pdb=s.complex_pdb,
                posebusters_valid=s.posebusters_valid,
                strain=s.strain,
                rmsd=s.rmsd,
            )
            for s in self._workflow.scores
        ]

    @property
    def conformers(self) -> list[str]:
        """UUIDs of optimized conformers."""
        return list(self._workflow.conformers)


class DoubleEndedTSSearchResult(WorkflowResult):
    """Result from a double-ended transition state search workflow."""

    _stjames_class = stjames.DoubleEndedTSSearchWorkflow

    @property
    def ts_guess_calculation_uuid(self) -> str | None:
        """UUID of the transition state guess calculation."""
        return self._workflow.ts_guess_calculation_uuid

    @property
    def distances(self) -> list[float] | None:
        """Path distances from reactant to product."""
        return list(self._workflow.distances) if self._workflow.distances else None


class ElectronicPropertiesResult(WorkflowResult):
    """Result from an electronic properties workflow."""

    _stjames_class = stjames.ElectronicPropertiesWorkflow


class FukuiResult(WorkflowResult):
    """Result from a Fukui index workflow."""

    _stjames_class = stjames.FukuiIndexWorkflow

    @property
    def global_electrophilicity_index(self) -> float | None:
        """Global electrophilicity index."""
        return self._workflow.global_electrophilicity_index

    @property
    def fukui_positive(self) -> list[float] | None:
        """Fukui f+ indices (electrophilic attack susceptibility)."""
        return list(self._workflow.fukui_positive) if self._workflow.fukui_positive else None

    @property
    def fukui_negative(self) -> list[float] | None:
        """Fukui f- indices (nucleophilic attack susceptibility)."""
        return list(self._workflow.fukui_negative) if self._workflow.fukui_negative else None

    @property
    def fukui_zero(self) -> list[float] | None:
        """Fukui f0 indices (radical attack susceptibility)."""
        return list(self._workflow.fukui_zero) if self._workflow.fukui_zero else None


class HydrogenBondBasicityResult(WorkflowResult):
    """Result from a hydrogen bond basicity workflow."""

    _stjames_class = stjames.HydrogenBondBasicityWorkflow


class IonMobilityResult(WorkflowResult):
    """Result from an ion mobility workflow."""

    _stjames_class = stjames.IonMobilityWorkflow

    @property
    def average_ccs(self) -> float | None:
        """Average collision cross section (Angstrom^2)."""
        return self._workflow.average_ccs

    @property
    def average_ccs_stdev(self) -> float | None:
        """Uncertainty in average CCS."""
        return self._workflow.average_ccs_stdev

    @property
    def conformer_ccs(self) -> list[float]:
        """Collision cross section per conformer (Angstrom^2)."""
        return list(self._workflow.conformer_ccs)

    @property
    def boltzmann_weights(self) -> list[float]:
        """Boltzmann weights for conformers."""
        return list(self._workflow.boltzmann_weights)


class IRCResult(WorkflowResult):
    """Result from an Intrinsic Reaction Coordinate (IRC) workflow."""

    _stjames_class = stjames.IRCWorkflow

    @property
    def starting_ts(self) -> str | None:
        """UUID of optimized TS before IRC."""
        return self._workflow.starting_TS

    @property
    def irc_forward(self) -> list[str]:
        """Forward IRC path UUIDs."""
        return list(self._workflow.irc_forward)

    @property
    def irc_backward(self) -> list[str]:
        """Backward IRC path UUIDs."""
        return list(self._workflow.irc_backward)


class MacropKaResult(WorkflowResult):
    """Result from a macropKa workflow."""

    _stjames_class = stjames.MacropKaWorkflow

    @property
    def isoelectric_point(self) -> float | None:
        """Isoelectric point (pH units)."""
        return self._workflow.isoelectric_point

    @property
    def solvation_energy(self) -> float | None:
        """Solvation energy (kcal/mol)."""
        return self._workflow.solvation_energy

    @property
    def kpuu_probability(self) -> float | None:
        """Probability that Kpuu >= 0.3."""
        return self._workflow.kpuu_probability

    @property
    def microstates(self) -> list[MacropKaMicrostate]:
        """All microstates."""
        return [
            MacropKaMicrostate(
                smiles=m.smiles,
                energy=m.energy,
                charge=m.charge,
            )
            for m in self._workflow.microstates
        ]

    @property
    def pka_values(self) -> list[MacropKaValue]:
        """Macroscopic pKa values."""
        return [
            MacropKaValue(
                initial_charge=v.initial_charge,
                final_charge=v.final_charge,
                pka=v.pKa,
            )
            for v in self._workflow.pKa_values
        ]

    @property
    def logd_by_ph(self) -> list[tuple[float, float]]:
        """Distribution constant by pH as (pH, logD) pairs."""
        return list(self._workflow.logD_by_pH)

    @property
    def aqueous_solubility_by_ph(self) -> list[tuple[float, float]]:
        """Aqueous solubility by pH as (pH, log(S)/L) pairs."""
        return list(self._workflow.aqueous_solubility_by_pH)


class MembranePermeabilityResult(WorkflowResult):
    """Result from a membrane permeability workflow."""

    _stjames_class = stjames.MembranePermeabilityWorkflow


class MSAResult(WorkflowResult):
    """Result from a Multiple Sequence Alignment (MSA) workflow."""

    _stjames_class = stjames.MSAWorkflow


class MultiStageOptResult(WorkflowResult):
    """Result from a multistage optimization workflow."""

    _stjames_class = stjames.MultiStageOptWorkflow


class NMRResult(WorkflowResult):
    """Result from a Nuclear Magnetic Resonance (NMR) workflow."""

    _stjames_class = stjames.NMRSpectroscopyWorkflow

    @property
    def chemical_shifts(self) -> list[float | None]:
        """Per-atom NMR chemical shifts (ensemble average)."""
        return list(self._workflow.chemical_shifts)

    @property
    def boltzmann_weights(self) -> list[float]:
        """Boltzmann weights for conformers."""
        return list(self._workflow.boltzmann_weights)

    @property
    def predicted_peaks(self) -> dict[int, list[NMRPeak]]:
        """Predicted NMR peaks by nucleus atomic number."""
        return {
            nucleus: [
                NMRPeak(
                    nucleus=p.nucleus,
                    shift=p.shift,
                    atom_indices=tuple(p.atom_indices),
                )
                for p in peaks
            ]
            for nucleus, peaks in self._workflow.predicted_peaks.items()
        }

    @property
    def symmetry_equivalent_nuclei(self) -> list[list[int]]:
        """0-indexed atoms which are equivalent to one another."""
        return self._workflow.symmetry_equivalent_nuclei


class pKaResult(WorkflowResult):
    """Result from a pKa workflow."""

    _stjames_class = stjames.pKaWorkflow

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


class PoseAnalysisMDResult(WorkflowResult):
    """Result from a Pose-Analysis Molecular Dynamics (MD) workflow."""

    _stjames_class = stjames.PoseAnalysisMolecularDynamicsWorkflow


class ProteinBinderDesignResult(WorkflowResult):
    """Result from a protein binder design workflow."""

    _stjames_class = stjames.ProteinBinderDesignWorkflow

    @property
    def generated_binders(self) -> list[ProteinBinder]:
        """Generated protein binder designs."""
        return [
            ProteinBinder(
                sequence=b.sequence,
                plddt=getattr(b, "plddt", None),
                ptm=getattr(b, "ptm", None),
                iptm=getattr(b, "iptm", None),
            )
            for b in self._workflow.generated_binders
        ]


class ProteinCofoldingResult(WorkflowResult):
    """Result from a protein cofolding workflow."""

    _stjames_class = stjames.ProteinCofoldingWorkflow

    @property
    def cofolding_results(self) -> list[CofoldingResult]:
        """Cofolding results."""
        return [
            CofoldingResult(
                plddt=getattr(r, "plddt", None),
                ptm=getattr(r, "ptm", None),
                iptm=getattr(r, "iptm", None),
            )
            for r in self._workflow.cofolding_results
        ]


class ProteinMDResult(WorkflowResult):
    """Result from a Protein Molecular Dynamics (MD) workflow."""

    _stjames_class = stjames.ProteinMolecularDynamicsWorkflow


class RedoxPotentialResult(WorkflowResult):
    """Result from a redox potential workflow."""

    _stjames_class = stjames.RedoxPotentialWorkflow

    @property
    def oxidation_potential(self) -> float | None:
        """Oxidation potential in V."""
        return self._workflow.oxidation_potential

    @property
    def reduction_potential(self) -> float | None:
        """Reduction potential in V."""
        return self._workflow.reduction_potential


class ScanResult(WorkflowResult):
    """Result from a scan workflow."""

    _stjames_class = stjames.ScanWorkflow

    @property
    def scan_points(self) -> list[str | None]:
        """UUIDs of scan point calculations."""
        return list(self._workflow.scan_points)


class SolubilityResult(WorkflowResult):
    """Result from an aqueous solubility workflow."""

    _stjames_class = stjames.SolubilityWorkflow

    @property
    def solubilities(self) -> list[SolubilityEntry]:
        """Solubility results per solvent."""
        return [
            SolubilityEntry(
                solvent=solvent,
                solubilities=tuple(result.solubilities),
                uncertainties=tuple(result.uncertainties),
            )
            for solvent, result in self._workflow.solubilities.items()
        ]

    @property
    def temperatures(self) -> list[float]:
        """Temperatures in Kelvin."""
        return list(self._workflow.temperatures)


class SpinStatesResult(WorkflowResult):
    """Result from a spin states workflow."""

    _stjames_class = stjames.SpinStatesWorkflow

    @property
    def spin_states(self) -> list[SpinState]:
        """List of spin states with energies."""
        return [
            SpinState(
                multiplicity=ss.multiplicity,
                energy=ss.energy,
            )
            for ss in self._workflow.spin_states
        ]

    @property
    def energies(self) -> list[float]:
        """Convenience: extract energies from spin_states."""
        return self._workflow.energies


class StrainResult(WorkflowResult):
    """Result from a strain workflow."""

    _stjames_class = stjames.StrainWorkflow

    @property
    def strain(self) -> float | None:
        """Computed strain energy (kcal/mol)."""
        return self._workflow.strain

    @property
    def conformers(self) -> list[str | None]:
        """UUIDs of conformers."""
        return list(self._workflow.conformers)


class TautomerResult(WorkflowResult):
    """Result from a tautomer search workflow."""

    _stjames_class = stjames.TautomerWorkflow

    @property
    def tautomers(self) -> list[Tautomer]:
        """List of tautomer structures."""
        return [
            Tautomer(
                energy=t.energy,
                weight=t.weight,
                predicted_relative_energy=t.predicted_relative_energy,
            )
            for t in self._workflow.tautomers
        ]


RESULT_REGISTRY: dict[str, type[WorkflowResult]] = {
    "admet": ADMETResult,
    "analogue_docking": AnalogueDockingResult,
    "basic_calculation": BasicCalculationResult,
    "batch_docking": BatchDockingResult,
    "bde": BDEResult,
    "conformer_search": ConformerSearchResult,
    "descriptors": DescriptorsResult,
    "docking": DockingResult,
    "double_ended_ts_search": DoubleEndedTSSearchResult,
    "electronic_properties": ElectronicPropertiesResult,
    "fukui": FukuiResult,
    "hydrogen_bond_basicity": HydrogenBondBasicityResult,
    "ion_mobility": IonMobilityResult,
    "irc": IRCResult,
    "macropka": MacropKaResult,
    "membrane_permeability": MembranePermeabilityResult,
    "msa": MSAResult,
    "multistage_opt": MultiStageOptResult,
    "nmr": NMRResult,
    "pka": pKaResult,
    "pose_analysis_md": PoseAnalysisMDResult,
    "protein_binder_design": ProteinBinderDesignResult,
    "protein_cofolding": ProteinCofoldingResult,
    "protein_md": ProteinMDResult,
    "redox_potential": RedoxPotentialResult,
    "scan": ScanResult,
    "solubility": SolubilityResult,
    "spin_states": SpinStatesResult,
    "strain": StrainResult,
    "tautomers": TautomerResult,
}


def create_result(
    workflow_data: dict,
    workflow_type: str,
    workflow_metadata: dict | None = None,
) -> WorkflowResult:
    """
    Factory function to create the appropriate result type for a workflow.

    :param workflow_data: The raw data dict from the workflow
    :param workflow_type: The workflow type string
    :param workflow_metadata: Optional metadata dict (elapsed, credits_charged, etc.)
    :return: A typed WorkflowResult subclass, or base WorkflowResult if unknown
    """
    result_class = RESULT_REGISTRY.get(workflow_type, WorkflowResult)
    return result_class(
        workflow_data=workflow_data,
        workflow_type=workflow_type,
        workflow_metadata=workflow_metadata or {},
    )
