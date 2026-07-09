import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

# BTK, catalytic Cys481. 4YHF already has a covalently-bound small-molecule inhibitor
# (residue 4C9) linked to Cys481 SG via the ligand's C1 atom.
protein = rowan.create_protein_from_pdb_id("4YHF")
protein = protein.select_chains(["A"])
protein.prepare(remove_heterogens=False)

# 0-based indices of Cys481 SG and the ligand's C1 atom, found by downloading the
# prepared PDB (protein.download_pdb_file()).
protein_reactive_atom_index = 1571
ligand_reactive_atom_index = 4492

workflow = rowan.submit_covalent_inhibitor_scan_workflow(
    protein=protein,
    protein_reactive_atom_index=protein_reactive_atom_index,
    ligand_reactive_atom_index=ligand_reactive_atom_index,
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
