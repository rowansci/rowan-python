"""Workflow submission and result types."""

# Import workflow-specific modules (each registers its result type)
from .admet import ADMETResult, submit_admet_workflow
from .analogue_docking import (
    AnalogueDockingResult,
    submit_analogue_docking_workflow,
)
from .base import (
    RESULT_REGISTRY,
    Workflow,
    WorkflowError,
    WorkflowResult,
    batch_submit_workflow,
    create_result,
    register_result,
    retrieve_workflow,
    submit_workflow,
)
from .basic_calculation import BasicCalculationResult, submit_basic_calculation_workflow
from .batch_docking import BatchDockingResult, submit_batch_docking_workflow
from .bde import BDEEntry, BDEResult, submit_bde_workflow
from .conformer_search import ConformerSearchResult, submit_conformer_search_workflow
from .descriptors import DescriptorsResult, submit_descriptors_workflow
from .docking import DockingResult, DockingScore, submit_docking_workflow
from .double_ended_ts_search import (
    DoubleEndedTSSearchResult,
    submit_double_ended_ts_search_workflow,
)
from .electronic_properties import (
    ElectronicPropertiesResult,
    submit_electronic_properties_workflow,
)
from .fukui import FukuiResult, submit_fukui_workflow
from .hydrogen_bond_basicity import (
    HydrogenBondBasicityResult,
    submit_hydrogen_bond_basicity_workflow,
)
from .ion_mobility import IonMobilityResult, submit_ion_mobility_workflow
from .irc import IRCResult, submit_irc_workflow
from .macropka import (
    MacropKaMicrostate,
    MacropKaResult,
    MacropKaValue,
    submit_macropka_workflow,
)
from .membrane_permeability import (
    MembranePermeabilityResult,
    submit_membrane_permeability_workflow,
)
from .msa import MSAResult, submit_msa_workflow
from .multistage_optimization import (
    MultiStageOptResult,
    submit_multistage_optimization_workflow,
)
from .nmr import NMRPeak, NMRResult, submit_nmr_workflow
from .pka import pKaMicrostate, pKaResult, submit_pka_workflow
from .pose_analysis_md import PoseAnalysisMDResult, submit_pose_analysis_md_workflow
from .protein_binder_design import (
    ProteinBinder,
    ProteinBinderDesignResult,
    submit_protein_binder_design_workflow,
)
from .protein_cofolding import (
    CofoldingResult,
    ProteinCofoldingResult,
    submit_protein_cofolding_workflow,
)
from .protein_md import ProteinMDResult, submit_protein_md_workflow
from .redox_potential import RedoxPotentialResult, submit_redox_potential_workflow
from .scan import ScanResult, submit_scan_workflow
from .solubility import SolubilityEntry, SolubilityResult, submit_solubility_workflow
from .spin_states import SpinState, SpinStatesResult, submit_spin_states_workflow
from .strain import StrainResult, submit_strain_workflow
from .tautomer_search import Tautomer, TautomerResult, submit_tautomer_search_workflow

# Backwards compatibility alias
AnalogueDockingScore = DockingScore

__all__ = [
    "RESULT_REGISTRY",
    "ADMETResult",
    "AnalogueDockingResult",
    "AnalogueDockingScore",
    "BDEEntry",
    "BDEResult",
    "BasicCalculationResult",
    "BatchDockingResult",
    "CofoldingResult",
    "ConformerSearchResult",
    "DescriptorsResult",
    "DockingResult",
    "DockingScore",
    "DoubleEndedTSSearchResult",
    "ElectronicPropertiesResult",
    "FukuiResult",
    "HydrogenBondBasicityResult",
    "IRCResult",
    "IonMobilityResult",
    "MSAResult",
    "MacropKaMicrostate",
    "MacropKaResult",
    "MacropKaValue",
    "MembranePermeabilityResult",
    "MultiStageOptResult",
    "NMRPeak",
    "NMRResult",
    "PoseAnalysisMDResult",
    "ProteinBinder",
    "ProteinBinderDesignResult",
    "ProteinCofoldingResult",
    "ProteinMDResult",
    "RedoxPotentialResult",
    "ScanResult",
    "SolubilityEntry",
    "SolubilityResult",
    "SpinState",
    "SpinStatesResult",
    "StrainResult",
    "Tautomer",
    "TautomerResult",
    "Workflow",
    "WorkflowError",
    "WorkflowResult",
    "batch_submit_workflow",
    "create_result",
    "pKaMicrostate",
    "pKaResult",
    "register_result",
    "retrieve_workflow",
    "submit_admet_workflow",
    "submit_analogue_docking_workflow",
    "submit_basic_calculation_workflow",
    "submit_batch_docking_workflow",
    "submit_bde_workflow",
    "submit_conformer_search_workflow",
    "submit_descriptors_workflow",
    "submit_docking_workflow",
    "submit_double_ended_ts_search_workflow",
    "submit_electronic_properties_workflow",
    "submit_fukui_workflow",
    "submit_hydrogen_bond_basicity_workflow",
    "submit_ion_mobility_workflow",
    "submit_irc_workflow",
    "submit_macropka_workflow",
    "submit_membrane_permeability_workflow",
    "submit_msa_workflow",
    "submit_multistage_optimization_workflow",
    "submit_nmr_workflow",
    "submit_pka_workflow",
    "submit_pose_analysis_md_workflow",
    "submit_protein_binder_design_workflow",
    "submit_protein_cofolding_workflow",
    "submit_protein_md_workflow",
    "submit_redox_potential_workflow",
    "submit_scan_workflow",
    "submit_solubility_workflow",
    "submit_spin_states_workflow",
    "submit_strain_workflow",
    "submit_tautomer_search_workflow",
    "submit_workflow",
]
