# Membrane permeability

## Input

A molecule. The accepted input type depends on the `method`:

- `gnn-mtl` requires a SMILES string passed as `initial_molecule`. Passing anything else raises `ValueError`.
- `pypermm` requires a 3D structure with coordinates (a `rowan.Molecule`, `stjames.Molecule`, or RDKit `Mol` carrying coordinates), not a SMILES string. A SMILES string, or a structure without coordinates (such as a bare `Chem.MolFromSmiles(...)`), raises `ValueError`.

## Example

```python
import rowan

folder = rowan.get_folder("examples")

smiles = "CC1=C(N=CN1)CSCCNC(=NC)NC#N"

# gnn-mtl: graph neural network, predicts experimental Caco-2 apparent permeability
gnn_wf = rowan.submit_membrane_permeability_workflow(
    smiles,
    method="gnn-mtl",
    folder=folder,
)

# pypermm: physics-based free-energy method, predicts intrinsic permeability for five membranes
pypermm_wf = rowan.submit_membrane_permeability_workflow(
    rowan.Molecule.from_smiles(smiles),   # pypermm needs a 3D structure with coordinates
    method="pypermm",
    folder=folder,
)

gnn_result = gnn_wf.result()
pypermm_result = pypermm_wf.result()

print(gnn_result.caco2_p_app)      # log10 Caco-2 apparent permeability P_app, in cm/s
print(pypermm_result.caco2_log_p)  # Caco-2 intrinsic permeability coefficient logP
```

## Methods

- `gnn-mtl` (default): a graph neural network (multi-task learning) that predicts the experimental Caco-2 apparent permeability `P_app`. Populates `caco2_p_app`.
- `pypermm`: PyPermm, a physics-based free-energy method that predicts intrinsic permeability coefficients (logP) for five membranes: plasma membrane, blood-brain barrier, Caco-2, black lipid membrane, and PAMPA. Also returns the free-energy profile across the membrane.

## Result fields

- `caco2_p_app`: log10 of the Caco-2 apparent permeability `P_app` in cm/s (gnn-mtl).
- `caco2_log_p`: Caco-2 intrinsic permeability coefficient logP (pypermm).
- `plasma_log_p`: plasma-membrane intrinsic permeability coefficient logP (pypermm).
- `bbb_log_p`: blood-brain-barrier intrinsic permeability coefficient logP (pypermm).
- `blm_log_p`: black lipid (bilayer) membrane intrinsic permeability coefficient logP (pypermm).
- `pampa_log_p`: PAMPA intrinsic permeability coefficient logP (pypermm).
- `energy_profile`: `(position (A), energy (kcal/mol))` pairs describing how the molecule's energy changes as it crosses the membrane (pypermm).
