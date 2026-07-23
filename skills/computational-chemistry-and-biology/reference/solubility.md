# Solubility

## Input

A SMILES string passed as `initial_smiles`. These models are 2D/SMILES-based: pass a SMILES string, not a 3D structure. To start from a `rowan.Molecule` or RDKit `Mol`, extract its SMILES first (`molecule.smiles` for a `rowan.Molecule`) and pass that string.

## Example

```python
import rowan

folder = rowan.get_folder("examples")

wf = rowan.submit_solubility_workflow(
    initial_smiles="C1CCOC(=O)C1=C[C@@H](OC(CC)CC)[C@H](NC(C)=O)[C@@H]([NH3+])C1CCC1",  # oseltamivir
    method="fastsolv",  # temperature-dependent, any solvent
    solvents=["ethanol", "toluene"],
    temperatures=[298.15, 323.15],
    folder=folder,
)

result = wf.result()
for entry in result.solubilities:
    for v in entry.values:
        print(
            f"  {entry.solvent} @ {v.temperature} K: {v.solubility:.2f} +/- {v.uncertainty} log(mol/L)"
        )
```

## Settings

- `method` (default `"fastsolv"`): solubility model. Choose by solvent: use `fastsolv` for any non-aqueous solvent, and for water use `kingfisher` (the recommended aqueous model) or `esol` (a simpler, lower-cost regression baseline).
  - `fastsolv`: ML-based solid solubility. Supports arbitrary solvents and temperatures — the go-to for non-aqueous solvents. Trained across roughly −30 to 130 °C, so keep temperatures in that range. Predicts water poorly, so don't use it for aqueous solubility.
  - `kingfisher`: ML-based aqueous solubility. Water only, 298.15 K only. The recommended model for aqueous solubility.
  - `esol`: ESOL regression for aqueous solubility. Water only, 298.15 K only. A simple, well-established baseline for water.
- `solvents` (default depends on method): list of solvent names or SMILES. Common names like `"ethanol"`, `"water"`, and `"thf"` are recognized, and for `fastsolv` any solvent SMILES is accepted. For `kingfisher` and `esol` it must be `["water"]` or `["O"]`. When omitted, `fastsolv` defaults to hexane, toluene, THF, ethyl acetate, ethanol, and acetonitrile, while `kingfisher` and `esol` default to water.
- `temperatures` (default depends on method): list of temperatures in Kelvin. Only `fastsolv` supports a temperature range — any temperatures are allowed (default `[273.15, 298.15, 323.15, 348.15, 373.15]`). `kingfisher` and `esol` are room-temperature only and must be `[298.15]`.

Solubility values are reported in log(mol/L).

## Result fields

- `solubilities`: one `SolubilityEntry` per solvent, each with `.solvent` and `.values`.
- each value (`.values`) carries `.temperature` (K), `.solubility` (log(mol/L)), and `.uncertainty` — the predicted standard deviation (populated for the ML models; `None` when a method does not report one).
