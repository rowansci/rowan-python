"""
Run a double-ended transition state search using Rowan.

See documentation at: https://docs.rowansci.com/science/workflows/double-ended-ts-search
"""

from stjames import Method, Molecule, Settings
from stjames.optimization.freezing_string_method import (
    FSMInterpolation,
    FSMOptimizationCoordinates,
    FSMSettings,
)

import rowan

# rowan.api_key = ""
HCN = Molecule.from_xyz(
    """\
H 0 0 -1.1
C 0 0 0
N 0 0 1.2""",
)
CNH = Molecule.from_xyz(
    """\
H 0 0 2.3
C 0 0 0
N 0 0 1.2""",
)
fsm_settings = FSMSettings(
    optimization_coordinates=FSMOptimizationCoordinates.CARTESIAN,
    interpolation_method=FSMInterpolation.REDUNDANT_INTERNAL_COORDINATES,
    min_num_nodes=7,
    num_interpolation_points=5,
    max_optimizer_iterations=3,
    max_line_search_steps=2,
    max_displacement=0.1,
)


# Run calculation remotely
result = rowan.submit_double_ended_ts_search_workflow(
    reactant=HCN,
    product=CNH,
    calculation_settings=Settings(method=Method.GFN2_XTB),
    search_settings=fsm_settings,
    optimize_inputs=True,
    optimize_ts=True,
    name="H-C≡N Isomerization",
)

print(f"View workflow privately at: https://labs.rowansci.com/workflow/{result.uuid}")
result.wait_for_result().fetch_latest(in_place=True)

assert result and result.data
print(result.data["forward_string_distances"])
print(result.data["backward_string_distances"])

print(result)
