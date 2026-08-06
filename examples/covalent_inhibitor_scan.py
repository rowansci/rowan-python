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

# Locate Cys481 SG and the 4C9 ligand's C1 atom in the prepared structure, then convert
# their atom serials to zero-based PDB record indices.
model = protein.data["models"][0]
cys481 = model["polymer"]["A"]["residues"]["A.481"]
ligand_4c9 = next(residue for residue in model["non_polymer"].values() if residue["name"] == "4C9")
protein_atom_serial = next(
    int(serial) for serial, atom in cys481["atoms"].items() if atom["name"] == "SG"
)
ligand_atom_serial = next(
    int(serial) for serial, atom in ligand_4c9["atoms"].items() if atom["name"] == "C1"
)
entities = [
    residue for chain in model["polymer"].values() for residue in chain["residues"].values()
]
for section in ("non_polymer", "water", "branched"):
    entities.extend(model.get(section, {}).values())
atom_serials = sorted(int(serial) for entity in entities for serial in entity.get("atoms", {}))
protein_reactive_atom_index = atom_serials.index(protein_atom_serial)
ligand_reactive_atom_index = atom_serials.index(ligand_atom_serial)

workflow = rowan.submit_covalent_inhibitor_scan_workflow(
    protein=protein.uuid,
    protein_reactive_atom_index=protein_reactive_atom_index,
    ligand_reactive_atom_index=ligand_reactive_atom_index,
    ligand_smiles=ligand_smiles,
    settings=rowan.CovalentInhibitorScanSettings(scan_num=4),
    name="BTK covalent inhibitor scan (Cys481, 4YHF)",
    folder=folder,
)

print(
    f"View workflow privately at: https://labs.rowansci.com/covalent-inhibitor-scan/{workflow.uuid}"
)

result = workflow.result()
print(result)

for distance, energy in result.get_energies():
    print(f"  distance={distance:.3f} Å  energy={energy}")
