"""Cofolding with a structural template.

Templates bias the prediction toward a known structure - useful when a
homolog or a previously solved conformation should guide the fold. They are
supported by Boltz-2 and OpenFold-3.

A template's `protein` is the UUID of a PDB record stored in Rowan. The
optional `max_distance` (in A) turns the template into a soft potential
bounding deviation from it. `max_distance` is Boltz-2 only, and when set it
requires `use_potentials=True`. There is no default - a template has no
distance bound unless you set one.

This example pulls a real structure (bovine trypsin, PDB 3PTB) into Rowan and
uses it as a template, submitting the same fold two ways:

1. OpenFold-3 with a plain template (no max_distance, no potentials).
2. Boltz-2 with a distance-bounded template (max_distance + use_potentials).
"""

import rowan

# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

# Bovine trypsin sequence, matching PDB 3PTB.
SEQUENCE = (
    "IVGGYTCGANTVPYQVSLNSGYHFCGGSLINSQWVVSAAHCYKSGIQVRLGEDNINVVEGNEQFISASKSIVHPSYNSNTLNN"
    "DIMLIKLKSAAYTSYDVPLGTQCLISGWGNTKSSGTSYPDVLKCLKAPILSDSSCKSAYPGQITSNMFCAGYLEGGKDSCQGD"
    "SGGPVVCSGKLQGIVSWGSGCAQKNKPGVYTKVCNYVSWIKQTIASN"
)

# Pull a real structure into Rowan and reference it by UUID.
template_protein = rowan.create_protein_from_pdb_id(name="3PTB trypsin template", code="3PTB")

# 1. OpenFold-3 with a plain template.
openfold_workflow = rowan.submit_protein_cofolding_workflow(
    initial_protein_sequences=[SEQUENCE],
    model=rowan.CofoldingModel.OPENFOLD_3,
    templates=[rowan.CofoldingTemplate(protein=template_protein.uuid)],
    num_samples=3,
    name="Trypsin (templated, OpenFold-3)",
    folder=folder,
)

# 2. Boltz-2 with a distance-bounded template. max_distance requires use_potentials.
boltz_workflow = rowan.submit_protein_cofolding_workflow(
    initial_protein_sequences=[SEQUENCE],
    model=rowan.CofoldingModel.BOLTZ_2,
    templates=[rowan.CofoldingTemplate(protein=template_protein.uuid, max_distance=2.0)],
    use_potentials=True,
    num_samples=3,
    name="Trypsin (templated, Boltz-2)",
    folder=folder,
)

for label, workflow in [("OpenFold-3", openfold_workflow), ("Boltz-2", boltz_workflow)]:
    print(f"\n{label}")
    print(f"  View at: https://labs.rowansci.com/protein-cofolding/{workflow.uuid}")
    for i, pred in enumerate(workflow.result().predictions):
        ptm = pred.scores.ptm if pred.scores else None
        print(f"  sample {i}: ptm={ptm}")
