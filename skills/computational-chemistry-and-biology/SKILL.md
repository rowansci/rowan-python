---
name: computational-chemistry-and-biology
description: Run computational chemistry and structural biology calculations on the Rowan platform through Rowan MCP tools or the `rowan` Python package. Use when the user wants to run Rowan workflows, or mentions Rowan, computational chemistry, structural biology, or molecular or protein simulation.
---

# Computational Chemistry and Biology with Rowan

Run simulations on the [Rowan](https://rowansci.com) platform through Rowan MCP tools or the `rowan` Python package. This skill provides shared scientific strategy and workflow guidance for both interfaces.

**Interfaces are for syntax; skill is for strategy.** Use this skill and its workflow references to choose scientifically appropriate methods and settings. Trust `discover_workflow` or the installed package when an exact parameter differs from a reference.

## Choose an interface

Both interfaces are first-class. Use the shared scientific guidance and workflow references below with either one.

- **MCP-connected chat:** when Rowan tools are available, read [MCP execution](reference/mcp_execution.md) before calling them. Do not install Python or request an API key; authentication belongs to the MCP connection. Use `discover_workflow` as the exact parameter contract.
- **Python code, scripts, or notebooks:** when the user requests Python or no Rowan MCP tools are available, read [Python SDK](reference/python_sdk.md) before writing or running code. Trust the installed package when its exact API differs from an example.

For project or folder management, also read [folders and projects](reference/folders_and_projects.md). For callbacks, read [webhooks](reference/webhooks.md).

## Molecule inputs

Structure-based workflows take a `rowan.Molecule`, not a bare SMILES string. Build one from a SMILES with `rowan.Molecule.from_smiles(...)` (RDKit embeds a 3D structure) or from coordinates with `rowan.Molecule.from_xyz(...)` / `from_xyz_file(...)`. A plain SMILES string is accepted only by the SMILES-based (2D) workflows. Each reference file states which case its workflow is:

- **SMILES (2D)**: predicts from the graph, ignores 3D, so a plain SMILES string works (ADMET, solubility, macropKa, the SMILES pKa/permeability methods).
- **3D from a SMILES is fine**: pass `rowan.Molecule.from_smiles(...)`; the embedded structure is built and optimized as needed (basic calculation, NMR, descriptors, redox, BDE, conformer search, and most QM workflows).
- **Real coordinates needed**: pass a structure with meaningful geometry via `rowan.Molecule.from_xyz(...)`, a prior result's `.molecule`, or a docked/crystal pose (electronic properties, interaction-energy decomposition, strain, IRC, double-ended TS, 3D pKa methods).

For higher rigor, start a 3D workflow from the conformer search workflow's lowest-energy conformer.

## Protein inputs

After loading a protein with `rowan.create_protein_from_pdb_id` or `rowan.upload_protein`, prepare it before any docking, MD, or FEP workflow. The default `protein.prepare()` fixes nonstandard residues, adds missing atoms, and adds hydrogens. See [reference/protein_prep.md](reference/protein_prep.md) for `prepare` and `validate_protein_forcefield`, options, and the common validation-failure fix.

## Workflows

Pick the workflow matching the task by reading the descriptions, then open its reference file (linked at the end of each entry) for recommended settings and an example. Python's universal submit arguments live in [Python SDK](reference/python_sdk.md) and aren't repeated per workflow; MCP parameters come from `discover_workflow`.

Not every account can run every workflow. MCP users must consult the `mcp_supported_workflows` field returned by `account_status`; Python users must consult `rowan.whoami().enabled_workflows` (a few workflows, such as FEP, are gated by plan).

- **ADMET**: predict ADME-Tox properties (absorption, distribution, metabolism, excretion, toxicity) with an ML model, for fast early developability, PK, and tox triage. See [reference/admet.md](reference/admet.md).
- **Analogue docking**: pose analogues of an already-bound reference ligand into consistent, analogous poses, to align and prepare a congeneric series, for example ahead of an RBFE screen. See [reference/analogue_docking.md](reference/analogue_docking.md).
- **Basic calculation**: compute core QM results for a single molecular or periodic (crystal/surface) structure, covering single-point energies, geometry and transition-state optimizations (including in solvent), and frequencies and thermochemistry, with DFT, NNP, or semiempirical methods. See [reference/basic_calculation.md](reference/basic_calculation.md).
- **Batch docking**: fast docking and scoring of many protein-ligand poses (AutoDock Vina or QVina2), keeping only the final scores rather than poses, for high-throughput virtual screening of large libraries. See [reference/batch_docking.md](reference/batch_docking.md).
- **Binding affinity**: SQM-based binding affinity for one or more protein–ligand poses — useful for post-docking triage before RBFE, or when ligand topologies are too dissimilar for FEP. See [reference/binding_affinity.md](reference/binding_affinity.md).
- **BDE**: homolytic bond-dissociation energies for selected bonds, to rank radical stability, find the weakest bonds, or flag likely metabolic or degradation soft spots. See [reference/bde.md](reference/bde.md).
- **Conformer search**: explore conformational space with fast methods, then rank the generated conformers with a more accurate final method. See [reference/conformer_search.md](reference/conformer_search.md).
- **Covalent inhibitor scan**: scan the distance between a reactive protein atom and a reactive ligand atom (covalent to non-covalent, or the reverse), to map out the reaction energetics of covalent bond formation. See [reference/covalent_inhibitor_scan.md](reference/covalent_inhibitor_scan.md).
- **Descriptors**: compute electronic and cheminformatic molecular descriptors, to featurize molecules for cheminformatics and machine learning. See [reference/descriptors.md](reference/descriptors.md).
- **Docking**: dock a ligand into a protein site and return the pose, refined with a constrained NNP optimization, plus optional conformer search (for a per-pose strain estimate) and pre-docking geometry optimization, for accurate single-ligand pose prediction. Standard and covalent docking are both supported. See [reference/docking.md](reference/docking.md).
- **Double-ended TS search**: find the transition state connecting a reactant and a product structure, for when you have both endpoints rather than a single TS guess. See [reference/double_ended_ts_search.md](reference/double_ended_ts_search.md).
- **Electronic properties**: compute molecular orbitals, electron density, electrostatic potential, atom-centered charges, bond orders, and multipole moments for a given structure (does not optimize geometry), to analyze electronic structure, reactivity, and charge distribution. See [reference/electronic_properties.md](reference/electronic_properties.md).
- **Fukui**: atom-centered Fukui indices (per-site reactivity toward nucleophiles, electrophiles, and radicals) plus the global electrophilicity index, to predict which sites react and how electrophilic a molecule is overall. See [reference/fukui.md](reference/fukui.md).
- **Hydrogen-bond acceptor/donor strength**: predict hydrogen-bond acceptor (pKBHX) and donor (pKα) strengths using neural network potentials and r2SCAN-3c DFT. See [reference/hydrogen_bond_donor_acceptor_strength.md](reference/hydrogen_bond_donor_acceptor_strength.md).
- **Interaction energy decomposition**: break the non-covalent interaction energy between two molecules into SAPT0 components (electrostatics, exchange, induction, dispersion), to quantify how much each contributes. See [reference/interaction_energy_decomposition.md](reference/interaction_energy_decomposition.md).
- **Ion mobility**: predict a molecule's rotationally averaged collision cross-section (CCS) in nitrogen via explicit trajectory-based scattering, accounting for conformers and protonation states. See [reference/ion_mobility.md](reference/ion_mobility.md).
- **IRC**: trace the steepest-descent reaction path in both directions from an optimized transition state to the reactant and product it connects. See [reference/irc.md](reference/irc.md).
- **macropKa**: enumerate a molecule's protonation microstates and conformers to predict its macroscopic pKa, microstate populations, and pH-dependent properties (isoelectric point, logD, aqueous solubility, blood-brain-barrier permeability). See [reference/macropka.md](reference/macropka.md).
- **Membrane permeability**: predict small-molecule permeability, either experimental Caco-2 apparent permeability via a graph neural network, or intrinsic permeability coefficients for five membranes (plasma, blood-brain barrier, Caco-2, black lipid membrane, PAMPA) via a physics-based free-energy method. See [reference/membrane_permeability.md](reference/membrane_permeability.md).
- **MSA**: generate and format multiple-sequence-alignment data for downstream use with co-folding models. See [reference/msa.md](reference/msa.md).
- **Multistage optimization**: quickly optimize a structure through successively more accurate levels of theory, yielding an accurate geometry and energy. See [reference/multistage_optimization.md](reference/multistage_optimization.md).
- **NMR**: predict a molecule's NMR spectrum. See [reference/nmr.md](reference/nmr.md).
- **pKa**: predict microscopic (per-site) pKa values, with a choice of ML-based and semiempirical methods. See [reference/pka.md](reference/pka.md).
- **Pocket detection**: detect potential binding pockets on a protein structure. See [reference/pocket_detection.md](reference/pocket_detection.md).
- **Pose-analysis MD**: run MD on a docked protein-ligand complex to test pose stability. See [reference/pose_analysis_md.md](reference/pose_analysis_md.md).
- **Protein binder design**: generate, filter, and rank protein binders for a given protein or small-molecule target. See [reference/protein_binder_design.md](reference/protein_binder_design.md).
- **Protein co-folding**: predict 3D structures of biomolecules and protein-ligand complexes using AlphaFold 3-style models (Boltz-2, Chai-1r, Boltz-1). See [reference/protein_cofolding.md](reference/protein_cofolding.md).
- **Protein MD**: run a molecular dynamics simulation on a protein. See [reference/protein_md.md](reference/protein_md.md).
- **RBFE graph**: take a list of ligands and generate the graph used by the relative binding free energy perturbation (FEP) workflow. Its edges are the ligand pairs FEP transforms between. See [reference/rbfe_graph.md](reference/rbfe_graph.md).
- **Redox potential**: predict a molecule's reduction or oxidation potential. See [reference/redox_potential.md](reference/redox_potential.md).
- **Relative binding free energy perturbation**: predict relative binding affinities across a congeneric ligand series using free energy perturbation (FEP), to rank analogs and prioritize which to synthesize. See [reference/relative_binding_free_energy_perturbation.md](reference/relative_binding_free_energy_perturbation.md).
- **Scan**: compute how a system's energy changes along a scan coordinate (bond, angle, or dihedral). See [reference/scan.md](reference/scan.md).
- **Solubility**: predict temperature-dependent solubility in a variety of solvents, plus room-temperature aqueous solubility from dedicated aqueous-solubility models. See [reference/solubility.md](reference/solubility.md).
- **Solvent-dependent conformers**: generate, optimize, and score conformers across multiple solvents to see how conformational preferences shift between them. See [reference/solvent_dependent_conformers.md](reference/solvent_dependent_conformers.md).
- **Spin states**: predict a molecule's lowest-energy spin state by optimizing it across different spin states. See [reference/spin_states.md](reference/spin_states.md).
- **Strain**: compute the strain energy of a specific pose by comparing it to the global-minimum conformer. See [reference/strain.md](reference/strain.md).
- **Tautomer search**: enumerate tautomers and rank their relative stability using neural network potentials or semiempirical methods. See [reference/tautomer_search.md](reference/tautomer_search.md).
