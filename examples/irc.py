from stjames import Molecule

import rowan

# Set ROWAN_API_KEY environment variable to your API key or set rowan.api_key directly
# rowan.api_key = "rowan-sk..."

workflow = rowan.submit_irc_workflow(
    initial_molecule=Molecule.from_xyz_lines(
        """7
SMILES `N=C([O-])[OH2+]`
N    -0.15519741  -1.36979175  -0.20679433
C     1.11565384  -1.23943631  -0.14797646
O     2.17614993  -1.72950370  -0.04017850
H    -0.55869366  -2.29559315  -0.23834737
O     1.02571386   0.42871733  -0.27925360
H    -0.09029954  -0.04166676  -0.31495768
H     1.26740151   0.88347299   0.53620841
""".splitlines()
    ),
    name="HNCO + H₂O - IRC",
    preopt=False,
)

print(f"View workflow privately at: https://labs.rowansci.com/workflow/{workflow.uuid}")
workflow.wait_for_result().fetch_latest(in_place=True)
print(workflow)
