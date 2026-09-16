# Basic Calculation

## Configuration

Use the exported enums instead of raw strings when setting the method and engine:

```python
import rowan

molecule = rowan.Molecule.from_smiles("CCO")
workflow = rowan.submit_basic_calculation_workflow(
    initial_molecule=molecule,
    tasks=["optimize"],
    method=rowan.Method.GFN2_XTB,
    engine=rowan.Engine.XTB,
)
```

Pass an [`OptimizationSettings`](settings.md) object or
an equivalent dictionary to `opt_settings`. `optimize_cell` defaults to `False`; set it to `True`
only when optimizing a periodic cell.

```python
opt_settings = rowan.OptimizationSettings(max_steps=200, optimize_cell=False)
```

::: rowan.workflows.basic_calculation
    handler: python
    options:
      show_source: false
      show_root_heading: false
      show_root_toc_entry: false
      members_order: source
      group_by_category: true
      filters: ["!^_"]
