"""Tests for scan energy and coordinate alignment."""

from unittest.mock import MagicMock, call

from pytest import MonkeyPatch, approx, mark

import rowan
import stjames
from rowan.workflows import scan as scan_module
from rowan.workflows.constants import HARTREE_TO_KCAL

MOLECULE = rowan.Molecule.from_xyz("2\nH2\nH 0 0 0\nH 0 0 0.74")
SETTINGS = rowan.ScanSettings(type="bond", atoms=[1, 2], start=0.7, stop=0.9, num=3)


@mark.parametrize(
    ("second_dimension", "repeats", "coordinates"),
    [([], 1, [0.7, 0.8, 0.9]), ([SETTINGS], 3, list(range(9)))],
)
@mark.parametrize(
    ("point_uuids", "expected_energies", "relative_energies"),
    [
        (["first", None, "last"], [-1, None, -0.99], [0, None, 0.01 * HARTREE_TO_KCAL]),
        ([None, None, None], [None, None, None], [None, None, None]),
        (["first", "unfinished", "last"], [-1, None, -0.99], [0, None, 0.01 * HARTREE_TO_KCAL]),
        ([None, "first", None], [None, -1, None], [None, 0, None]),
    ],
)
def test_scan_energies_preserve_missing_point_positions(
    monkeypatch: MonkeyPatch,
    second_dimension: list[rowan.ScanSettings],
    repeats: int,
    coordinates: list[float],
    point_uuids: list[str | None],
    expected_energies: list[float | None],
    relative_energies: list[float | None],
) -> None:
    """Keep missing energies aligned with scan coordinates and flattened grid indices."""
    workflow = stjames.ScanWorkflow(
        initial_molecule=MOLECULE.to_stjames(),
        scan_settings=SETTINGS,
        scan_settings_2d=second_dimension,
        calc_settings=rowan.Settings(method="gfn2_xtb"),
        scan_points=point_uuids * repeats,
    )
    calculations = {
        uuid: rowan.Calculation(
            uuid=uuid,
            molecules=[
                rowan.Molecule.from_stjames(
                    MOLECULE.to_stjames().model_copy(update={"energy": energy})
                )
            ],
        )
        for uuid, energy in {"first": -1, "unfinished": None, "last": -0.99}.items()
    }
    retrieve = MagicMock(side_effect=calculations.__getitem__)
    monkeypatch.setattr(scan_module, "retrieve_calculation", retrieve)
    result = scan_module.ScanResult(
        workflow_data=workflow.model_dump(mode="json"),
        workflow_type="scan",
        workflow_uuid="workflow",
        complete=False,
    )

    energies = result.get_energies()
    relative = result.get_energies(relative=True)

    assert [coordinate for coordinate, _ in energies] == approx(coordinates)
    assert [coordinate for coordinate, _ in relative] == approx(coordinates)
    assert [energy for _, energy in energies] == approx(expected_energies * repeats)
    assert [energy for _, energy in relative] == approx(relative_energies * repeats)
    assert retrieve.call_args_list == [call(uuid) for uuid in point_uuids * repeats if uuid]
