# LogP

## Input

A SMILES string passed as `initial_smiles`. This workflow is SMILES-based: pass a SMILES string, not a 3D structure. To start from a `rowan.Molecule` or RDKit `Mol`, extract its SMILES first (`molecule.smiles` for a `rowan.Molecule`) and pass that string. The `cosmors` method generates its own conformers internally, so it also takes only a SMILES.

## Example

```python
import rowan

folder = rowan.get_folder("examples")

wf = rowan.submit_logp_workflow(
    initial_smiles="C1CCOC(=O)C1=C[C@@H](OC(CC)CC)[C@H](NC(C)=O)[C@@H]([NH3+])C1CCC1",  # oseltamivir
    method="chemprop_sangster2026",
    folder=folder,
)

result = wf.result()
print(result.logp)
```

## Settings

- `method` (default `"chemprop_sangster2026"`): logP prediction method.
  - `chemprop_sangster2026`: chemprop v2 D-MPNN trained on experimental octanol/water logP from the Sangster dataset. Fast, and the recommended default for drug-like molecules.
  - `crippen`: RDKit's Wildman-Crippen atom-contribution model. Essentially free, so useful for filtering large libraries, but it is an additive scheme that misses conformational and intramolecular effects.
  - `cosmors`: physics-based. Optimizes a conformer ensemble (GFN2-xTB/ALPB, then g-xTB/CPCM-X in water), runs COSMO-RS on the Boltzmann-significant conformers, and returns the Boltzmann-weighted logP. Much slower than the other two, but it uses no training data, so reach for it on chemistry poorly represented in experimental logP datasets.

## Result fields

- `logp`: predicted base-10 octanol/water partition coefficient.

logP describes the neutral species. For the pH-dependent equivalent (logD) on an ionizable compound, use the macropKa workflow.
