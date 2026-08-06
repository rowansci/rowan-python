import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

# TG2, catalytic Cys277, acrylamide inhibitor (k_inact record P21980-00).
tg2_inhibitor = rowan.Molecule.from_smiles(
    "C=CC(=O)NCCCC[C@H](NC(=O)Cc1ccc(Cl)c(Cl)c1)C(=O)N1CCN(C(=O)c2cccc3ccccc23)CC1"
)

protein = rowan.create_protein_from_pdb_id("2Q3Z")
# Chain A is TG2; chain X is a small covalently-bound peptide in the crystal, not
# part of the protein. `protein.chains` order isn't guaranteed, so select by name.
protein = protein.select_chains(["A"])
preparation_workflow = rowan.submit_protein_preparation_workflow(
    protein=protein.uuid,
    name="Prepare TG2",
    folder=folder,
)
protein = preparation_workflow.result().get_prepared_protein()

# Pocket is [[center_x, center_y, center_z], [size_x, size_y, size_z]] in Å.
center = [-1.079, -3.081, 18.122]
size = [22.22, 14.08, 21.74]

# 0-based, all-atom (including hydrogens) indices of the reacting atoms. The ligand's
# acrylamide terminal =CH2 carbon is atom 0. Locate Cys277 SG from the prepared structure
# so the protein index remains correct if preparation changes atom ordering.
covalent_ligand_atom_index = 0
model = protein.data["models"][0]
cys277 = model["polymer"]["A"]["residues"]["A.277"]
cys277_sg_serial = next(
    int(serial) for serial, atom in cys277["atoms"].items() if atom["name"] == "SG"
)
entities = [
    residue for chain in model["polymer"].values() for residue in chain["residues"].values()
]
for section in ("non_polymer", "water", "branched"):
    entities.extend(model.get(section, {}).values())
atom_serials = sorted(int(serial) for entity in entities for serial in entity.get("atoms", {}))
covalent_protein_atom_index = atom_serials.index(cys277_sg_serial)
gnina_settings = rowan.GninaSettings(
    scoring_function="gnina_cnn",
    covalent_ligand_atom_index=covalent_ligand_atom_index,
    covalent_protein_atom_index=covalent_protein_atom_index,
)

workflow = rowan.submit_docking_workflow(
    protein.uuid,
    pocket=[center, size],
    initial_molecule=tg2_inhibitor,
    docking_settings=gnina_settings,
    name="TG2 covalent docking (Cys277, 2Q3Z)",
    folder=folder,
)

print(f"View workflow privately at: https://labs.rowansci.com/docking/{workflow.uuid}")

result = workflow.result()
print(result)

for i, score in enumerate(result.scores):
    print(f"  Pose {i}: score={score.score:.3f}  posebusters_valid={score.posebusters_valid}")
