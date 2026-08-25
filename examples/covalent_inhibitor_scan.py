import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

# BTK, catalytic Cys481. 4YHF already has a covalently-bound small-molecule inhibitor
# (residue 4C9) linked to Cys481 SG via the ligand's C1 atom.
ligand_smiles = "CC(C)(C)C[C@@H](C#N)C(=O)N1CCC[C@H](C1)n2nc(c3ccc(Oc4ccccc4)cc3)c5c(N)ncnc25"

protein = rowan.create_protein_from_pdb_id("4YHF")
protein = protein.select_chains(["A"])
preparation_workflow = rowan.submit_protein_preparation_workflow(
    protein=protein.uuid,
    retain_non_polymer={"4C9": ligand_smiles},
    name="Prepare BTK inhibitor complex",
    folder=folder,
)
protein = preparation_workflow.result().get_prepared_protein()

# Protein preparation normalizes 4YHF's residue numbering: Cys481 becomes residue 101,
# while the retained 4C9 ligand remains residue 701.
protein_reactive_atom_index = protein.get_atom_index(chain="A", residue=101, atom="SG")
ligand_reactive_atom_index = protein.get_atom_index(
    chain="A", residue=701, atom="C1", entity_type="non_polymer"
)

workflow = rowan.submit_covalent_inhibitor_scan_workflow(
    protein=protein.uuid,
    protein_reactive_atom_index=protein_reactive_atom_index,
    ligand_reactive_atom_index=ligand_reactive_atom_index,
    reactant_smiles=ligand_smiles,
    name="BTK covalent inhibitor scan (Cys481, 4YHF)",
    folder=folder,
)

print(
    f"View workflow privately at: https://labs.rowansci.com/covalent-inhibitor-scan/{workflow.uuid}"
)

result = workflow.result()
print(result)

for distance, free_energy in result.get_energies():
    print(f"  distance={distance:.3f} Å  free_energy={free_energy} kcal/mol")
