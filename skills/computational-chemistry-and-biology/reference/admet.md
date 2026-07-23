# ADMET

## Input

A SMILES string passed as `initial_smiles`. These models are 2D/SMILES-based: pass a SMILES string, not a 3D structure. To start from a `rowan.Molecule` or RDKit `Mol`, extract its SMILES first (`molecule.smiles` for a `rowan.Molecule`) and pass that string.

## Example

```python
import rowan

folder = rowan.get_folder("examples")

oseltamivir = "C1CCOC(=O)C1=C[C@@H](OC(CC)CC)[C@H](NC(C)=O)[C@@H]([NH3+])C1CCC1"

wf = rowan.submit_admet_workflow(
    initial_smiles=oseltamivir,
    folder=folder,
)

result = wf.result()
print(result.properties)  # dict of property name -> predicted value
```

## Settings

This workflow takes no scientific settings beyond the input molecule. The result exposes a single field:

- `properties`: a dict of predicted properties spanning physicochemical, absorption, distribution, metabolism, excretion, and toxicity endpoints.

Global models of ADME/toxicity are not substitutes for experimental measurements. Predicting macroscopic properties like solubility, permeability, and toxicity is incredibly difficult, and while today's models are good enough to give decent guesses most of the time, they often still get things wrong. Verify any crucial predictions with experimental data.
