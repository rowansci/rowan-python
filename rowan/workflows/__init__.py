"""Workflow submission and result types.

Submit workflows and retrieve typed results. See examples/ directory
for 45+ detailed examples of all workflow types.

Basic pattern:
    result = workflows.submit_<workflow_type>(molecule).result()
"""

from stjames import (
    ETKDGSettings,
    LyrebirdSettings,
    MonteCarloMultipleMinimumSettings,
    SolventModel,
    iMTDGCSettings,
)

from .admet import ADMETResult, submit_admet_workflow
from .analogue_docking import (
    AnalogueDockingResult,
    submit_analogue_docking_workflow,
)
from .base import (
    RESULT_REGISTRY,
    Solvent,
    Workflow,
    WorkflowError,
    WorkflowResult,
    batch_poll_status,
    batch_submit_workflow,
    create_result,
    list_workflows,
    register_result,
    retrieve_calculation_molecules,
    retrieve_workflow,
    retrieve_workflows,
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
from .hydrogen_bond_donor_acceptor_strength import (
    HydrogenBondAcceptorSite,
    HydrogenBondBasicityResult,
    HydrogenBondDonorAcceptorStrengthResult,
    HydrogenBondDonorSite,
    submit_hydrogen_bond_basicity_workflow,
    submit_hydrogen_bond_donor_acceptor_strength_workflow,
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
