import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples/phenol-pka")

phenols_to_compute = {
    "p-nitrophenol": "Oc1ccc(N(=O)=O)cc1",
    "p-(trifluoromethyl)phenol": "Oc1ccc(C(F)(F)F)cc1",
    "methyl p-hydroxybenzoate": "COC(=O)c1ccc(O)cc1",
    "p-fluorophenol": "Oc1ccc(F)cc1",
    "p-chlorophenol": "Oc1ccc(Cl)cc1",
    "phenol": "Oc1ccccc1",
    "p-cresol": "Oc1ccc(C)cc1",
    "p-methoxyphenol": "COc1ccc(O)cc1",
    "p-(dimethylamino)phenol": "CN(C)c1ccc(O)cc1",
}

workflows = []
for name, smiles in phenols_to_compute.items():
    workflows.append(
        rowan.submit_pka_workflow(
            rowan.Molecule.from_smiles(smiles),
            name=f"pKa {name}",
            folder=folder,
            deprotonate_elements=[8],
        )
    )


workflow_results = [w.result() for w in workflows]

print(
    [(w.name, r.conjugate_bases[0].pka) for w, r in zip(workflows, workflow_results, strict=True)]
)
