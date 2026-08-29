# Batch Docking

Use `num_poses_to_save` to retain the best pose for that many top-scoring compounds, and set
`run_mmgbsa=True` to refine those saved poses. `result.refined_scores` maps every input SMILES
to a `DockingScore` for a retained compound or `None` otherwise.

::: rowan.workflows.batch_docking
    handler: python
    options:
      show_source: false
      show_root_heading: false
      show_root_toc_entry: false
      members_order: source
      group_by_category: true
      filters: ["!^_"]
