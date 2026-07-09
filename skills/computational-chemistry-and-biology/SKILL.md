---
name: computational-chemistry-and-biology
description: Run computational chemistry and structural biology calculations on the Rowan platform via the `rowan` Python package. Use when the user wants to run Rowan workflows, or mentions Rowan, computational chemistry, structural biology, or molecular or protein simulation.
---

# Computational Chemistry and Biology with Rowan

Run simulations on the [Rowan](https://rowansci.com) platform via the `rowan` Python package. Every workflow follows the same flow: build an input, submit a workflow, then retrieve the result.

**Package is for syntax, skill is for strategy.** For exact arguments, defaults, and valid values, introspect the installed package with `python3 -c "import rowan; help(rowan.submit_pka_workflow)"` (any name works). This skill tells you which workflow to use and recommends settings. If a value here is rejected, trust the package.

## Setup

Do this before calling the Rowan API. Use the resulting env's python for every later step.

**1. Install `rowan`** (PyPI: `rowan-python`). Skip if `python3 -c "import rowan"` already succeeds. Otherwise add it through the project's env manager so it survives syncs. A bare `pip install` into a managed env gets wiped:

- **uv** (`uv.lock`): `uv add rowan-python`, then run with `uv run python`
- **pixi** (`pixi.toml`): `pixi add --pypi rowan-python`, then `pixi run python`
- **poetry** (`poetry.lock`): `poetry add rowan-python`, then `poetry run python`
- **plain venv or system**: `python3 -m pip install rowan-python`

**2. Set the key and verify.** Put `ROWAN_API_KEY=your-key-here` in a project `.env` (gitignored, and create a key at https://labs.rowansci.com/account), then run the preflight. It finds the nearest `.env`, validates the key against the API, and prints remaining credits (workflows cost credits). Run it from your project root, where `.env` lives, and point it at this skill's base dir rather than `cd`-ing into the skill:

`python3 "<skill_dir>/scripts/check_env.py"`

It exits with actionable guidance on each failure: a **missing** key (add it via `echo 'ROWAN_API_KEY=...' >> .env`), a **rejected** key (a 401, meaning invalid or expired, so regenerate it and stop to tell the user), or an **unreachable** API. Since shell state doesn't persist between commands, load `.env` in the same command as every later call, like `[ -f .env ] && set -a && source .env && set +a; <python ...>`, or set `rowan.api_key = os.environ["ROWAN_API_KEY"]` in Python.

## Run & retrieve

Every workflow is submitted and retrieved the same way, regardless of type. `submit_*_workflow(...)` returns a `Workflow` immediately, and results are fetched separately.

```python
import rowan

folder = rowan.get_folder("my-project")          # path created if it doesn't exist
wf = rowan.submit_<workflow>_workflow(            # see its reference file (Workflows below) for scientific args
    ...,                                          # workflow-specific arguments
    folder=folder,
    name="my run",
)
result = wf.result()                              # blocks until done, returns a typed result
```

**Universal submit arguments.** On top of its scientific arguments, every `submit_*_workflow` accepts these, and they are not repeated in the per-workflow files:

- `folder=` (or `folder_uuid=`) sets where the run lands. `rowan.get_folder("a/b/c")` returns the folder, creating the path if missing. For a specific project or benchmark, make a dedicated folder (e.g. `get_folder("project-name")`, with subpaths like `"project-name/docking"` for sub-experiments) and route all its runs there, so results stay grouped and easy to retrieve later. To browse or navigate an existing folder tree, or switch projects, see [reference/folders_and_projects.md](reference/folders_and_projects.md).
- `name=` sets the label shown in the Rowan UI.
- `max_credits=` caps the credits (≈ minutes of compute) a run may use — its runtime limit. Hitting it, or running out of credits, stops the run **permanently** (methods like DFT can't resume), so set it high enough to finish.
- `webhook_url=` makes Rowan POST the result there on completion. See [reference/webhooks.md](reference/webhooks.md) if you're setting up webhooks.
- `is_draft=True` stages the run without starting it, for estimation (below).

**Three ways to retrieve results:**

- **Block (default):** `result = wf.result()` waits (polling every 5 s) and returns a typed result. It raises `rowan.WorkflowError` if the run **failed** or was **stopped**, so catch it to report failures.
- **Partial or progress:** `wf.result(wait=False)` returns whatever is ready now. `for r in wf.stream_result():` yields partial results each poll, then the final one.
- **Fire-and-forget:** keep `wf.uuid`, let the session end, and reconnect later with `rowan.retrieve_workflow(uuid).result()`. `rowan.list_workflows()` lists recent runs. The UUID from `labs.rowansci.com/calculation/<uuid>` is a workflow UUID — use `retrieve_workflow`, not `retrieve_calculation`.

**Inline calls vs. scripts.** Quick questions or retrieving a single result are fine inline with submit-and-wait (`.result()`). For long-running or multi-step work like experiments or benchmarks, put the calls in a script for reproducibility, and prefer a fire-and-forget run plus `retrieve_workflow(uuid)` later or a `webhook_url=` callback rather than blocking the session.

**Estimate cost and runtime before committing.** Submit as a draft (nothing runs), inspect, then commit or discard:

```python
draft = rowan.submit_<workflow>_workflow(..., folder=folder, is_draft=True)
print(draft.dispatch_info())   # DispatchInfo(to_be_dispatched, compute_hardware, estimated_runtime_minutes)
draft.submit_draft()           # commit to running it  (or draft.delete() to discard)
```

**Status and control** (non-blocking): `wf.done()` and `wf.get_status()` check state without waiting, `wf.stop()` cancels a running workflow, and `wf.delete()` removes it and its data.

**Resubmit and re-run.** Re-run from a settings dict with the generic `rowan.submit_workflow(workflow_type, workflow_data, initial_molecule=..., folder=folder)`. Pass a prior run's `result.data` unchanged for an identical re-run (insulated from later default or preset changes), edit the dict to change settings, or change `workflow_type` (with matching `workflow_data`) to submit as a different workflow.

**Resubmit with a perturbed geometry.** `mol.perturb()` adds small Gaussian noise to break symmetry or escape a stuck optimization. `mol.displace_along_mode(mode, displacement)` shifts atoms along a vibrational normal mode — for a transition state, pick the imaginary mode (negative frequency) to slide toward reactant or product; requires `frequencies=True` in a prior calculation. See `examples/resubmit_with_perturbations.py`.

## Molecule inputs

Structure-based workflows take a `rowan.Molecule`, not a bare SMILES string. Build one from a SMILES with `rowan.Molecule.from_smiles(...)` (RDKit embeds a 3D structure) or from coordinates with `rowan.Molecule.from_xyz(...)` / `from_xyz_file(...)`. A plain SMILES string is accepted only by the SMILES-based (2D) workflows. Each reference file states which case its workflow is:

- **SMILES (2D)**: predicts from the graph, ignores 3D, so a plain SMILES string works (ADMET, solubility, macropKa, the SMILES pKa/permeability methods).
- **3D from a SMILES is fine**: pass `rowan.Molecule.from_smiles(...)`; the embedded structure is built and optimized as needed (basic calculation, NMR, descriptors, redox, BDE, conformer search, and most QM workflows).
- **Real coordinates needed**: pass a structure with meaningful geometry via `rowan.Molecule.from_xyz(...)`, a prior result's `.molecule`, or a docked/crystal pose (electronic properties, interaction-energy decomposition, strain, IRC, double-ended TS, 3D pKa methods).

For higher rigor, start a 3D workflow from the conformer search workflow's lowest-energy conformer.

## Protein inputs

After loading a protein with `rowan.create_protein_from_pdb_id` or `rowan.upload_protein`, prepare it before any docking, MD, or FEP workflow. The default `protein.prepare()` fixes nonstandard residues, adds missing atoms, and adds hydrogens. See [reference/protein_prep.md](reference/protein_prep.md) for `prepare` and `validate_protein_forcefield`, options, and the common validation-failure fix.

## Workflows

Pick the workflow matching the task by reading the descriptions, then open its reference file (linked at the end of each entry) for recommended settings and an example. Universal submit arguments (`folder`, `name`, `max_credits`, and so on) live in **Run & retrieve** above and aren't repeated per workflow.

Not every key can run every workflow — `rowan.whoami().enabled_workflows` lists the types yours can submit (a few, e.g. FEP, are gated by plan).

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
