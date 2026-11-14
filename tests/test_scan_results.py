"""Tests for scan energy and coordinate alignment."""

from unittest.mock import MagicMock

import pytest
import stjames

import rowan
from rowan.workflows import scan as scan_module
from rowan.workflows.constants import HARTREE_TO_KCAL


@pytest.mark.parametrize("relative", [False, True])
@pytest.mark.parametrize("two_dimensional", [False, True])
def test_scan_energies_preserve_missing_point_positions(
    monkeypatch: pytest.MonkeyPatch, relative: bool, two_dimensional: bool
) -> None:
    """Keep missing calculations aligned with scan coordinates and flattened grid indices."""
    molecule = rowan.Molecule.from_xyz("2\n\nH 0 0 0\nH 0 0 0.74")
    settings = rowan.ScanSettings(type="bond", atoms=[1, 2], start=0.7, stop=0.9, num=3)
    point_uuids = [
        f"point-{index}" if index % 3 != 1 else None for index in range(9 if two_dimensional else 3)
    ]
    workflow = stjames.ScanWorkflow(
        initial_molecule=molecule.to_stjames(),
        scan_settings=settings,
        scan_settings_2d=settings if two_dimensional else [],
        calc_settings=rowan.Settings(method="gfn2_xtb"),
        scan_points=point_uuids,
    )
    calculations = {
        uuid: rowan.Calculation(
            uuid=uuid,
            molecules=[
                rowan.Molecule.from_stjames(
                    molecule.to_stjames().model_copy(
                        update={"energy": -1.0 if index % 3 == 0 else -0.99}
                    )
                )
            ],
        )
        for index, uuid in enumerate(point_uuids)
        if uuid is not None
    }
    retrieve = MagicMock(side_effect=calculations.__getitem__)
    monkeypatch.setattr(scan_module, "retrieve_calculation", retrieve)
    result = scan_module.ScanResult(
        workflow_data=workflow.model_dump(mode="json"),
        workflow_type="scan",
        workflow_uuid="workflow",
        complete=False,
    )

    energies = result.get_energies(relative=relative)

    coordinates = list(range(9)) if two_dimensional else [0.7, 0.8, 0.9]
    expected_energies = [0.0, None, 0.01 * HARTREE_TO_KCAL] if relative else [-1.0, None, -0.99]
    expected_energies *= 3 if two_dimensional else 1
    assert [coordinate for coordinate, _ in energies] == pytest.approx(coordinates)
    assert [energy for _, energy in energies] == pytest.approx(expected_energies)
    assert all(args.args[0] is not None for args in retrieve.call_args_list)
