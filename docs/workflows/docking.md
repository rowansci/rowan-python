# Docking

The docking workflow supports Vina docking and both noncovalent and covalent gnina docking.
Passing `GninaSettings` selects gnina; it does not by itself enable covalent docking.

## Covalent docking

Set both covalent atom indices on `GninaSettings` to form a bond between a known ligand atom and
protein atom. Covalent gnina docking requires `scoring_function="vina"`.

Prepare the protein first, then resolve the reactive protein atom from the prepared structure.
Protein preparation can change atom ordering and residue numbering.

Supply the ligand in its expected post-reaction, covalently bound topology; gnina does not infer the
reaction. For a Michael acceptor `C=CC(=O)NR`, use the hydrogen-capped product `CCC(=O)NR` and
select the terminal β-carbon as the covalent ligand atom.

```python
prepared_protein = preparation_workflow.result().get_prepared_protein()
reactive_protein_atom_index = prepared_protein.get_atom_index(
    chain="A", residue=reactive_residue, atom="SG"
)
settings = rowan.GninaSettings(
    scoring_function="vina",
    covalent_ligand_atom_index=reactive_ligand_atom_index,
    covalent_protein_atom_index=reactive_protein_atom_index,
)
workflow = rowan.submit_docking_workflow(
    prepared_protein.uuid,
    pocket=[center, size],
    initial_molecule=ligand,
    docking_settings=settings,
)
```

Both indices are zero-based all-atom indices, including hydrogens. See
`examples/covalent_docking.py` for a complete TG2 example.

PoseBusters validation is skipped for covalent poses. Their `posebusters_valid` value is `None`,
meaning not evaluated rather than failed; do not use it to reject covalent poses.

::: rowan.workflows.docking
    handler: python
    options:
      show_source: false
      show_root_heading: false
      show_root_toc_entry: false
      members_order: source
      group_by_category: true
      filters: ["!^_"]
