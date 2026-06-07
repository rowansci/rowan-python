"""
Rank a conformer ensemble you already have, using the Rowan API (screen-only mode).

Conformer search normally *generates* conformers and then optimizes, deduplicates,
and ranks them. If you already have 3D geometries -- from another tool (RDKit,
CREST, OpenBabel), a crystal or MD ensemble, an SDF on disk, or a previous
workflow -- you can skip generation by passing them as `initial_conformers`. Rowan
then runs only the screening half: optimize each at a consistent level of theory,
deduplicate, and rank by energy.

The conformers must be a genuine ensemble of one molecule with identical atom
ordering (they are compared atom-by-atom during deduplication). Reading them from a
single multi-conformer SDF, as below, guarantees that. Leave `conf_gen_settings` as
None.

See documentation at: https://docs.rowansci.com/science/workflows/conformers
"""

from pathlib import Path

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

# Read a multi-conformer SDF: eight conformers of ibuprofen, one molecule, consistent
# atom ordering across every record.
sdf_path = Path(__file__).parent / "data" / "ibuprofen_conformers.sdf"
conformers = rowan.Molecule.molecules_from_sdf(sdf_path)
print(f"read {len(conformers)} conformers from {sdf_path.name}")

workflow = rowan.submit_conformer_search_workflow(
    initial_conformers=conformers,
    conf_gen_settings=None,  # screen-only: skip generation, just optimize/dedup/rank
    name="Conformer Screen",
    folder=folder,
)

print(f"View workflow privately at: https://labs.rowansci.com/conformer-search/{workflow.uuid}")

result = workflow.result()
print(result)
# e.g. <ConformerSearchResult conformers=3>  (8 supplied -> deduplicated and ranked)

print("relative energies (kcal/mol):", [round(e, 2) for e in result.get_energies(relative=True)])
