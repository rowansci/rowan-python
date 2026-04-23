import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."

folder = rowan.get_folder("examples")

# Naphthalene-guanine non-covalent complex.
# Fragment 1 (atoms 1-18): naphthalene
# Fragment 2 (atoms 19-34): guanine
naphthalene_guanine_xyz = """\
34
name: Naphthalene-Guanine; charge: 0; multiplicity: 1;
C    -0.33944955   0.04784594   0.19447983
C    -0.35035696  -1.34375866   0.21601688
C    -1.56463757  -2.03222175   0.20796367
C    -2.78421442  -1.33642508   0.17822905
C    -4.01455643  -2.01360436   0.16975329
C    -5.21789589  -1.30676969   0.14010616
C    -5.20698562   0.08483354   0.11856891
C    -3.99270597   0.77329497   0.12662541
C    -2.77312954   0.07749639   0.15634444
C    -1.54278996   0.75467638   0.16482524
H     0.60532603   0.58494882   0.20070576
H     0.58589158  -1.89514251   0.23912042
H    -1.55294445  -3.11959000   0.22511859
H    -4.04330411  -3.10066600   0.18627183
H    -6.16267958  -1.84385641   0.13388886
H    -6.14323201   0.63621823   0.09544986
H    -4.00440588   1.86066393   0.10947873
H    -1.53740705   1.82973733   0.14323021
C    -4.78577410  -1.46541214   3.35535139
N    -3.77820183  -2.30853939   3.38755834
C    -2.64295229  -1.56778043   3.38787147
C    -2.99685089  -0.21618136   3.37082725
N    -4.36711799  -0.16139552   3.34009334
C    -1.98635390   0.79359056   3.30859618
O    -2.24061928   1.99196665   3.27852227
N    -0.71077365   0.26554727   3.26147078
C    -0.42452852  -1.09861474   3.17847860
N    -1.35369611  -2.02443841   3.27766533
N     0.87675342  -1.47664023   2.86846462
H    -5.82024848  -1.77873102   3.30060440
H    -4.91422115   0.69500666   3.22422884
H     0.03948138   0.93736360   3.18012383
H     1.57067595  -0.75990301   2.72540239
H     1.09557491  -2.30845030   2.31265810
"""

dimer = rowan.Molecule.from_xyz(naphthalene_guanine_xyz)

workflow = rowan.submit_interaction_energy_decomposition_workflow(
    initial_molecule=dimer,
    fragment1_indices=list(range(1, 19)),  # atoms 1-18: naphthalene
    folder=folder,
    name="Naphthalene-Guanine SAPT0",
)

print(f"View at: https://labs.rowansci.com/interaction-energy-decomposition/{workflow.uuid}")

result = workflow.result()
print(result)
# <InteractionEnergyDecompositionResult total_interaction_energy=-3.462>

print(f"Total:         {result.total_interaction_energy:>8.3f} kcal/mol")
# Total:            -3.462 kcal/mol
print(f"Electrostatic: {result.electrostatic_interaction_energy:>8.3f} kcal/mol")
# Electrostatic:   -15.450 kcal/mol
print(f"Exchange:      {result.exchange_interaction_energy:>8.3f} kcal/mol")
# Exchange:         38.280 kcal/mol
print(f"Dispersion:    {result.dispersion_interaction_energy:>8.3f} kcal/mol")
# Dispersion:      -23.695 kcal/mol
print(f"Induction:     {result.induction_interaction_energy:>8.3f} kcal/mol")
# Induction:        -2.597 kcal/mol
