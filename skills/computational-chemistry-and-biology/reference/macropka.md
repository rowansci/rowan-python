# macropKa

## Input

A SMILES string passed as `initial_smiles`. This is a 2D/SMILES-based workflow: it takes a SMILES string, not a 3D structure, and generates its own conformers. To start from a `rowan.Molecule` or RDKit `Mol`, extract its SMILES first (`molecule.smiles` for a `rowan.Molecule`) and pass that string. It enumerates the molecule's protonation microstates and conformers, then predicts macroscopic pKa values, microstate populations, the isoelectric point, pH-dependent logD and aqueous solubility, and (when `compute_solvation_energy=True`) a Kpuu-based blood-brain-barrier permeability probability.

## Example

```python
import rowan

folder = rowan.get_folder("examples")

wf = rowan.submit_macropka_workflow(
    # oseltamivir
    initial_smiles="C1CCOC(=O)C1=C[C@@H](OC(CC)CC)[C@H](NC(C)=O)[C@@H]([NH3+])C1CCC1",
    folder=folder,
)

result = wf.result()
print(result.pka_values)  # macroscopic pKa values; see Result fields for the rest
```

## Settings

- `min_pH` (default `0`) and `max_pH` (default `14`): pH range over which microstate weights, logD, and solubility are computed.
- `min_charge` (default `-2`) and `max_charge` (default `2`): range of net charges to enumerate microstates over. `min_charge` must be less than `max_charge`.
- `compute_aqueous_solubility` (default `True`): predict pH-dependent aqueous solubility. Non-ideal behavior such as aggregation is not modeled.
- `compute_solvation_energy` (default `False`): run a conformer search and compute the solvation energy, used to estimate Kpuu (unbound brain-to-plasma partition, a blood-brain-barrier permeability measure).

## Result fields

- `isoelectric_point`: the pH at which the molecule carries no net charge.
- `pka_values`: macroscopic pKa values, each as a charge transition with `initial_charge`, `final_charge`, and `pKa`.
- `microstates`: enumerated microstates, each with `smiles`, `energy` (free energy), and `charge`.
- `microstate_weights_by_ph`: `(pH, [weights])` pairs giving the population of each microstate at each pH.
- `logd_by_ph`: `(pH, logD)` pairs (water/octanol distribution coefficient).
- `aqueous_solubility_by_ph`: `(pH, log(S)/L)` pairs.
- `kpuu_probability`: probability that Kpuu is at least 0.3. Populated when `compute_solvation_energy=True`.
- `solvation_energy`: solvation energy in kcal/mol.

Macroscopic pKa values come from the Starling model (a retrained Uni-pKa model). On standard pKa benchmarks (SAMPL6/7/8 and the Novartis acid/base sets), its accuracy is comparable to state-of-the-art tools such as Uni-pKa, ChemAxon, and Epik, roughly 0.7 to 1.1 pKa units.
