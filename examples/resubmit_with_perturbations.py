"""
Resubmit a completed workflow with a perturbed structure.

Two strategies: random noise to break symmetry, or displacement along a
vibrational mode to follow a reaction coordinate or escape a stuck geometry.
"""

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

# --- Option 1: random noise ---
wf = rowan.retrieve_workflow("your-workflow-uuid")
mol = wf.result().molecule

perturbed_mol = mol.perturb()
resubmit = rowan.submit_multistage_optimization_workflow(
    initial_molecule=perturbed_mol,
    name="cyclobutane opt - perturbed resubmit",
    folder=folder,
)
print(f"https://labs.rowansci.com/multistage-opt/{resubmit.uuid}")

# --- Option 2: displace along a vibrational mode ---
# Requires a prior frequency calculation. Imaginary modes have negative frequency.
ts_wf = rowan.retrieve_workflow("your-ts-freq-workflow-uuid")
ts_mol = ts_wf.result().molecule

imaginary_mode = next(m for m in ts_mol.vibrational_modes if m.frequency < 0)

displaced_mol = ts_mol.displace_along_mode(mode=imaginary_mode, displacement=0.3)
resubmit = rowan.submit_multistage_optimization_workflow(
    initial_molecule=displaced_mol,
    name="TS resubmit - displaced along imaginary mode",
    folder=folder,
)
print(f"https://labs.rowansci.com/multistage-opt/{resubmit.uuid}")
