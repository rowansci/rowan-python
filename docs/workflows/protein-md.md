# Protein MD

Protein MD accepts selectable small-molecule, protein, and water force fields. Results expose
binder RMSD, MM/GBSA scores, mean-structure UUIDs, and medoid frame indices per trajectory. Use
`get_mean_structure`, `download_mean_structure`, or `download_medoid_structure` to retrieve
representative structures.

For protein structures containing non-polymer ligands, pass their SMILES through
`small_molecules` and select the residues to analyze through `rowan.Binder`.

Use `implicit_solvent="gbn2"` with `pressure_atm`, `nonbonded_cutoff`, and `water_buffer`
set to `None`. GBn2 uses nonperiodic dynamics and `ionic_strength_M` for salt screening.
REST2 and solvent saving are unsupported.

Pass `enhanced_sampling=rowan.REST2Settings()` to enable whole-solute REST2. The default
eight-state ladder runs from 298.15 to 450 K with exchange attempts every 1 ps. Custom
`effective_temperatures` must start at the workflow's `temperature`.
`num_trajectories` counts independent exchange ensembles, and `simulation_time_ns` applies
to each replica, so computational cost includes every state in the ladder. Each ensemble
returns one trajectory following the physical (λ=1) state through exchanges. Structural
analyses use these frames; their sequence cannot be used to infer kinetics or residence times.
Leave `enhanced_sampling=None` for conventional MD.

::: rowan.workflows.protein_md
    handler: python
    options:
      show_source: false
      show_root_heading: false
      show_root_toc_entry: false
      members_order: source
      group_by_category: true
      filters: ["!^_"]
