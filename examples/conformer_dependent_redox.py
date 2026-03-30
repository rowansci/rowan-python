import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples/conformer-dependent-redox")

workflow = rowan.submit_conformer_search_workflow(
    initial_molecule=rowan.Molecule.from_smiles("CC(C)Cc1ccc(C(=O)c2ccc(O)cc2)cc1"),
    folder=folder,
)
print(f"View workflow privately at: https://labs.rowansci.com/workflow/{workflow.uuid}")
csearch_result = workflow.result()

redox_potential_workflows = []

for conformer in csearch_result.conformer_uuids[:10]:
    uuid = conformer[0]
    molecule = rowan.retrieve_calculation_molecules(uuid)[0]
    rowan_molecule = rowan.Molecule.model_validate(molecule)
    redox_potential_workflows.append(
        rowan.submit_redox_potential_workflow(
            rowan_molecule,
            reduction=True,
            oxidation=True,
            folder=folder,
        )
    )

redox_results = [w.result() for w in redox_potential_workflows]

print([(r.oxidation_potential, r.reduction_potential) for r in redox_results])
