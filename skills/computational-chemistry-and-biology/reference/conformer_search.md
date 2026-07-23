# Conformer search

## Input

Two mutually exclusive input modes:

- **Generation** (`initial_molecule`): generate and rank conformers from a single structure. OpenConf (the default) and ETKDG accept a SMILES string or a 3D `rowan.Molecule`; MCMM and iMTDGC require a 3D `rowan.Molecule` and raise a `ValueError` on a SMILES.
- **Screening** (`initial_conformers`): rank a supplied ensemble without generation. Requires 3D `rowan.Molecule` structures (no SMILES) and `initial_molecule` left unset.

## Example

```python
import rowan

folder = rowan.get_folder("examples")

wf = rowan.submit_conformer_search_workflow(
    initial_molecule="CCOCC",  # SMILES works because the default generator is OpenConf
    final_method="aimnet2_wb97md3",
    folder=folder,
)

result = wf.result()
print(result.get_energies(relative=True))  # relative energies (kcal/mol), lowest first
best = result.get_conformers(1)[0]  # lowest-energy conformer, a rowan.Molecule
```

`get_conformers(n)` returns the `n` lowest-energy conformers as `rowan.Molecule` objects (omit `n` for all), and `get_energies(relative=True)` gives their relative energies in kcal/mol. Both are energy-ordered, so the lowest-energy conformer is index 0, ready to feed into a downstream 3D workflow. `get_conformers` makes one API call per conformer, so request only as many as you need.

## Settings

- `final_method` (default `aimnet2_wb97md3`): method used to rank the final conformers.
- `solvent` (default gas phase): implicit solvent; introspect `rowan.Solvent` for the available solvents.
- `transition_state` (default `False`): optimize each conformer toward a transition state rather than a minimum. Start from a fully optimized transition state, not a guess and not a SMILES. Conformer generation will otherwise relax the transition state to a local minimum, so pass constraints in `conf_gen_settings` to hold the reacting atoms in place during generation; the post-generation TS optimizations then run without constraints.
- `conf_gen_settings` (default none, currently resolves to OpenConf): override conformer generation. Pick and tune a generator by passing its settings object: `rowan.OpenConfSettings`, `rowan.ETKDGSettings`, `rowan.MonteCarloMultipleMinimumSettings` (MCMM), or `rowan.iMTDGCSettings` (iMTDGC). OpenConf (the default) is Rowan's conformer generator; ETKDG is a fast RDKit knowledge-based generator; iMTDGC runs CREST metadynamics for more thorough sampling of flexible molecules, at higher cost. Each has its own tuning fields; run `help()` on the one you want. `max_confs` is common to every generator and caps how many conformers are kept after generation.
- `multistage_opt_settings` (default none): override the full optimization-and-ranking stack as a `MultiStageOptSettings`. When supplied, it takes precedence over `final_method`, `solvent`, and `transition_state`, so those are ignored. Omit it and set `final_method`/`solvent`/`transition_state` instead unless you need a stack they cannot express.
- `conformer_clustering_settings` (default none): thin the generated ensemble to representative conformers before the expensive optimization, as a `rowan.ConformerClusteringSettings`. Clustering follows the ReSCoSS protocol (k-means on 3D-shape descriptors, keeping the lowest-energy conformers from each cluster). Raise `num_clusters` for more diversity across the output, raise `conformers_per_cluster` for more coverage within each region. Not supported with `initial_conformers` (there is no generated ensemble to cluster), which raises a `ValueError`.
- `initial_conformers` (default none): pre-generated 3D conformers to rank directly (screen-only mode: skip generation, just optimize, deduplicate, and rank). Mutually exclusive with `initial_molecule`, and leave `conf_gen_settings` at none. The conformers must be a genuine ensemble of one molecule with identical atom ordering (they are compared atom-by-atom during deduplication); load a multi-conformer SDF with `rowan.Molecule.molecules_from_sdf(path)` to guarantee that.
