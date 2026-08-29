# Protein MD

Protein MD accepts selectable small-molecule, protein, and water force fields. Results expose
binder RMSD, MM/GBSA scores, mean-structure UUIDs, and medoid frame indices per trajectory. Use
`get_mean_structure`, `download_mean_structure`, or `download_medoid_structure` to retrieve
representative structures.

For protein structures containing non-polymer ligands, pass their SMILES through
`small_molecules` and select the residues to analyze through `rowan.Binder`.

::: rowan.workflows.protein_md
    handler: python
    options:
      show_source: false
      show_root_heading: false
      show_root_toc_entry: false
      members_order: source
      group_by_category: true
      filters: ["!^_"]
