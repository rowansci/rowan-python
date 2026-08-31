"""Cofolding with modified polymers and a covalent bond constraint."""

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable.
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

protein = rowan.ProteinSequence(
    sequence="ASA",
    modifications=[rowan.ResidueModification(position=1, ccd="SEP")],
)
dna = rowan.DNASequence(
    sequence="AC",
    modifications=[rowan.NucleotideModification(position=1, ccd="5MC")],
)
rna = rowan.RNASequence(
    sequence="AU",
    modifications=[rowan.NucleotideModification(position=1, ccd="PSU")],
)

protein_atom = rowan.ConstraintTarget(
    input_type="protein", input_index=0, token_index=1, atom_name="OG"
)
ligand_atom = rowan.ConstraintTarget(input_type="ligand", input_index=0, token_index=0)
covalent_bond = rowan.BondConstraint(atom_1=protein_atom, atom_2=ligand_atom)

workflow = rowan.submit_protein_cofolding_workflow(
    initial_protein_sequences=[protein],
    initial_dna_sequences=[dna],
    initial_rna_sequences=[rna],
    initial_smiles_list=["CBr"],
    bond_constraints=[covalent_bond],
    model=rowan.CofoldingModel.BOLTZ_2,
    num_samples=1,
    name="Modified polymers with covalent constraint",
    folder=folder,
)

print(f"View workflow privately at: https://labs.rowansci.com/protein-cofolding/{workflow.uuid}")
result = workflow.result()
print(result)
for i, prediction in enumerate(result.predictions):
    print(f"  sample {i}: scores={prediction.scores}")
