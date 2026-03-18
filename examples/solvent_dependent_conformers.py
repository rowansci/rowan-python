import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."

folder = rowan.get_folder("examples")

alanine_dipeptide = rowan.Molecule.from_smiles("CC(=O)N[C@@H](C)C(=O)NC")

workflow = rowan.submit_solvent_dependent_conformers_workflow(
    initial_molecule=alanine_dipeptide,
    folder=folder,
    name="Alanine Dipeptide Solvent-Dependent Conformers",
)

print(f"View at: https://labs.rowansci.com/workflow/{workflow.uuid}")

result = workflow.result()
print(result)
# <SolventDependentConformersResult conformers=8 solvents=5>

# Transfer free energies (relative to lowest solvent)
for solvent, dg in result.relative_free_energy_by_solvent.items():
    print(f"  {solvent.value}: ΔG = {dg:.2f} kcal/mol")

# Per-solvent ensemble properties
for solvent, props in result.per_solvent_properties.items():
    print(
        f"  {solvent.value}: SASA={props.solvent_accessible_surface_area:.1f} Å², "
        f"r_g={props.radius_of_gyration:.2f} Å"
    )

# Per-conformer populations in water
for i, conf in enumerate(result.conformers):
    pop = conf.population_by_solvent.get(rowan.Solvent.WATER, 0)
    dg = conf.relative_free_energy_by_solvent.get(rowan.Solvent.WATER, float("nan"))
    print(f"  Conformer {i + 1}: ΔG={dg:.2f} kcal/mol, population={pop * 100:.1f}% (water)")
