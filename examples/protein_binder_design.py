# ruff: noqa: E501
"""
Design protein binders using BoltzGen via the Rowan API.

This workflow generates, filters, and ranks protein designs that bind
to given protein or small-molecule targets.

Available protocols:
- "protein-anything": Design proteins to bind proteins (default)
- "peptide-anything": Design peptides to bind proteins
- "protein-small_molecule": Design proteins to bind small molecules
- "nanobody-anything": Design single-domain antibodies

Sequence format for designable regions:
- "140..180": Design a protein of 140-180 residues
- "9A9C": Design 9 residues, fix an alanine, design 9 more, fix a cysteine
- "MKTAYIAKQ": Fixed amino-acid sequence (no designable region)

See documentation at: https://docs.rowansci.com/science/workflows/protein-binder-design
"""

import rowan

# Set ROWAN_API_KEY environment variable to your API key or set rowan.api_key directly
# rowan.api_key = "rowan-sk..."

# Design a protein to bind brilacidin (antimicrobial small molecule)
# Based on: https://github.com/HannesStark/boltzgen/blob/main/example/protein_binding_small_molecule/brilacidin.yaml
binder_design_input = {
    "protein_entities": [
        {
            "id": "A",
            "sequence": "140..180",  # design a 140-180 residue protein binder
        },
    ],
    "ligand_entities": [
        {
            "id": "B",
            # brilacidin - antimicrobial small molecule
            "smiles": "C1CNC[C@@H]1OC2=C(C=C(C=C2NC(=O)C3=CC(=NC=N3)C(=O)NC4=CC(=CC(=C4O[C@@H]5CCNC5)NC(=O)CCCCN=C(N)N)C(F)(F)F)C(F)(F)F)NC(=O)CCCCN=C(N)N",
        },
    ],
}

workflow = rowan.submit_protein_binder_design_workflow(
    binder_design_input=binder_design_input,
    protocol="protein-small_molecule",
    num_designs=10,  # number of designs to generate
    budget=2,  # number of designs to return after diversity filtering
    name="Brilacidin binder design",
)

print(
    f"View workflow privately at: https://labs.rowansci.com/protein-binder-design/{workflow.uuid}"
)
workflow.wait_for_result().fetch_latest(in_place=True)
print(workflow.data)
