# Types

Value types shared across many workflows — the vocabulary you use for `preset`, `method`, and
structural inputs regardless of which workflow you're submitting.

Legacy enum members retained for parsing historical workflows are omitted below.

Enum members and their accepted string values are listed below. Import these types from `rowan`,
not `stjames`.

## Methods

::: stjames.method.Method
    handler: python
    options:
      show_source: false
      show_bases: false
      show_root_heading: false
      show_root_toc_entry: false
      show_if_no_docstring: true
      separate_signature: false
      members_order: source
      group_by_category: true
      filters:
        - "!^_"
        - "!^default_engine$"
        - "!^(MACE_MP_0|MACE_MP_0B2_L|EGRET_1|EGRET_1E|EGRET_1T)$"
        - "!^SMIRNOFF_2_(0_0|2_1)_AMBER_AM1BCC$"

## Engines

::: stjames.engine.Engine
    handler: python
    options:
      show_source: false
      show_bases: false
      show_root_heading: false
      show_root_toc_entry: false
      show_if_no_docstring: true
      separate_signature: false
      members_order: source
      group_by_category: true
      filters:
        - "!^_"
        - "!^(EGRET|MACE|TERACHEM)$"

## Tasks

::: stjames.task.Task
    handler: python
    options:
      show_source: false
      show_bases: false
      show_root_heading: false
      show_root_toc_entry: false
      show_if_no_docstring: true
      separate_signature: false
      members_order: source
      group_by_category: true
      filters:
        - "!^_"
        - "!^STRESS$"

## Corrections

::: stjames.correction.Correction
    handler: python
    options:
      show_source: false
      show_bases: false
      show_root_heading: false
      show_root_toc_entry: false
      show_if_no_docstring: true
      separate_signature: false
      members_order: source
      group_by_category: true
      filters: ["!^_"]

## Modes

::: stjames.mode.Mode
    handler: python
    options:
      show_source: false
      show_bases: false
      show_root_heading: false
      show_root_toc_entry: false
      show_if_no_docstring: true
      separate_signature: false
      members_order: source
      group_by_category: true
      filters: ["!^_"]

## Solvent

::: stjames.solvent
    handler: python
    options:
      show_source: false
      show_root_heading: false
      show_root_toc_entry: false
      members_order: source
      group_by_category: true
      filters: ["!^_"]

## Constraints

::: stjames.constraint
    handler: python
    options:
      show_source: false
      show_root_heading: false
      show_root_toc_entry: false
      members_order: source
      group_by_category: true
      filters: ["!^_"]
