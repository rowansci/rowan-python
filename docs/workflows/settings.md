# Settings

::: stjames.settings
    handler: python
    options:
      show_source: false
      show_root_heading: false
      show_root_toc_entry: false
      members_order: source
      group_by_category: true
      filters: ["!^_"]
      members:
        - Settings

## Basis sets

Pass a basis set name as a string for ordinary calculations. Use `rowan.BasisSet` when configuring
atom- or element-specific overrides.

::: stjames.basis_set
    handler: python
    options:
      show_source: false
      show_root_heading: false
      show_root_toc_entry: false
      members_order: source
      group_by_category: true
      filters: ["!^_"]

::: stjames.opt_settings
    handler: python
    options:
      show_source: false
      show_root_heading: false
      show_root_toc_entry: false
      members_order: source
      group_by_category: true
      filters: ["!^_"]

::: stjames.pbc_dft_settings
    handler: python
    options:
      show_source: false
      show_root_heading: false
      show_root_toc_entry: false
      members_order: source
      group_by_category: true
      filters: ["!^_"]

::: stjames.excited_state_settings
    handler: python
    options:
      show_source: false
      show_root_heading: false
      show_root_toc_entry: false
      members_order: source
      group_by_category: true
      filters: ["!^_"]
