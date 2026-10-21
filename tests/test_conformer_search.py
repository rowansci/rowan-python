"""Tests for conformer-search workflow submission."""

from collections.abc import Iterator
from contextlib import contextmanager
from unittest.mock import MagicMock

import pytest
import stjames

import rowan
from rowan.workflows import conformer_search as conformer_search_module


def _workflow_response() -> dict[str, object]:
    """Build a minimal workflow API response."""
    return {
        "name": "test",
        "uuid": "workflow-uuid",
        "created_at": "2026-09-15T00:00:00Z",
        "object_status": 0,
        "parent_uuid": "folder-uuid",
        "notes": "",
        "starred": False,
        "public": False,
        "object_type": "conformer_search",
        "object_data": {},
        "email_when_complete": False,
        "credits_charged": 0,
        "object_logfile": "",
    }


def _mock_submission(monkeypatch: pytest.MonkeyPatch) -> MagicMock:
    """Replace the conformer-search API client and return its mock client."""
    client = MagicMock()
    client.post.return_value.json.return_value = _workflow_response()

    @contextmanager
    def mock_api_client() -> Iterator[MagicMock]:
        yield client

    monkeypatch.setattr(conformer_search_module, "api_client", mock_api_client)
    return client


def test_conformer_constraints_are_released_for_ts_refinement(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Keep generation constrained while releasing TS refinement constraints."""
    client = _mock_submission(monkeypatch)
    constraint = rowan.Constraint(constraint_type="freeze_atoms", atoms=[1, 2])
    conf_gen_settings = rowan.iMTDGCSettings(constraints=[constraint])
    optimization_stage = stjames.Settings(
        method=stjames.Method.OMOL25_CONSERVING_S,
        tasks=[stjames.Task.OPTIMIZE],
        opt_settings={
            "transition_state": True,
            "constraints": [constraint],
        },
    )
    multistage_opt_settings = stjames.MultiStageOptSettings(
        optimization_settings=[optimization_stage]
    )

    rowan.submit_conformer_search_workflow(
        initial_molecule=rowan.Molecule.from_xyz("2\n\nH 0 0 0\nH 0 0 0.74"),
        conf_gen_settings=conf_gen_settings,
        multistage_opt_settings=multistage_opt_settings,
    )

    workflow_data = client.post.call_args.kwargs["json"]["workflow_data"]
    assert workflow_data["conf_gen_settings"]["constraints"] == [
        {"constraint_type": "freeze_atoms", "atoms": [1, 2]}
    ]
    assert (
        workflow_data["multistage_opt_settings"]["optimization_settings"][0]["opt_settings"][
            "constraints"
        ]
        == []
    )
    assert optimization_stage.opt_settings.constraints == [constraint]


ENERGIES = [0.0, 0.001]


def _conformer_molecule(
    charge: float, excitation: float, oscillator: float, rotatory: float
) -> rowan.Molecule:
    """H2 conformer carrying charges, one excited state, and one IR mode."""
    return rowan.Molecule.from_stjames(
        stjames.Molecule(
            charge=0,
            multiplicity=1,
            atoms=[
                stjames.Atom(atomic_number=1, position=(0, 0, 0)),
                stjames.Atom(atomic_number=1, position=(0, 0, 0.74)),
            ],
            mulliken_charges=[charge, -charge],
            excited_state_energies=[-1.0, -1.0 + excitation],
            oscillator_strengths=[0.0, oscillator],
            rotatory_strengths=[0.0, rotatory],
            vibrational_modes=[
                stjames.VibrationalMode(
                    frequency=frequency,
                    reduced_mass=1.0,
                    force_constant=1.0,
                    displacements=[(0, 0, 1), (0, 0, -1)],
                    ir_intensity=100.0,
                )
                for frequency in (-30.0, 4400.0)
            ],
        )
    )


def _ensemble_result(
    monkeypatch: pytest.MonkeyPatch, ensemble: stjames.ConformerProperties
) -> conformer_search_module.ConformerSearchResult:
    """Two-conformer result whose calculations are served from memory."""
    molecules = {
        "calc-0": _conformer_molecule(0.2, 0.2, 0.4, 10.0),
        "calc-1": _conformer_molecule(0.6, 0.3, 0.8, -20.0),
    }
    monkeypatch.setattr(
        conformer_search_module,
        "retrieve_calculation",
        lambda uuid: rowan.Calculation(uuid=uuid, molecules=[molecules[uuid]]),
    )
    workflow = stjames.ConformerSearchWorkflow(
        initial_smiles="[HH]",
        conf_gen_settings=rowan.ETKDGSettings(),
        conformer_uuids=[["calc-0"], ["calc-1"]],
        energies=ENERGIES,
        conformer_properties=[
            rowan.ConformerProperties(solvent_accessible_surface_area=100, dipole=1),
            rowan.ConformerProperties(radius_of_gyration=2, dipole=3),
        ],
        ensemble_properties=ensemble,
    )
    return conformer_search_module.ConformerSearchResult(
        workflow_data=workflow.model_dump(mode="json", serialize_as_any=True),
        workflow_type="conformer_search",
        workflow_uuid="workflow-uuid",
    )


def test_ensemble_properties_at_298_k_are_rowans(monkeypatch: pytest.MonkeyPatch) -> None:
    """Return Rowan's ensemble properties at 298.15 K without fetching calculations."""
    ensemble = rowan.ConformerProperties(dipole=1.5, mulliken_charges=[0.1, -0.1])
    result = _ensemble_result(monkeypatch, ensemble)
    monkeypatch.setattr(
        conformer_search_module, "retrieve_calculation", MagicMock(side_effect=AssertionError)
    )

    assert result.ensemble_properties == ensemble
    assert result.get_ensemble_properties() == ensemble
    assert result.dipoles == [1, 3]


@pytest.mark.parametrize("with_ecd", [True, False])
def test_ensemble_properties_reweight_at_temperature(
    monkeypatch: pytest.MonkeyPatch, with_ecd: bool
) -> None:
    """Reweight scalars, charges, and stick spectra at another temperature."""
    ensemble = rowan.ConformerProperties(
        mulliken_charges=[0.1, -0.1],
        excitation_energies=[0.2],
        oscillator_strengths=[0.1],
        rotatory_strengths=[1.0] if with_ecd else None,
        vibrational_frequencies=[4400.0],
        ir_intensities=[100.0],
    )
    result = _ensemble_result(monkeypatch, ensemble)
    w0, w1 = result.get_boltzmann_weights(1000)
    assert w0 + w1 == pytest.approx(1)
    assert w0 < result.get_boltzmann_weights()[0]

    properties = result.get_ensemble_properties(1000)
    assert properties.solvent_accessible_surface_area == pytest.approx(100)
    assert properties.radius_of_gyration == pytest.approx(2)
    assert properties.dipole == pytest.approx((w0 * 1 + w1 * 9) ** 0.5, abs=1e-3)
    assert properties.mulliken_charges == pytest.approx(
        [w0 * 0.2 + w1 * 0.6, -(w0 * 0.2 + w1 * 0.6)], abs=1e-6
    )
    assert properties.excitation_energies == pytest.approx([0.2, 0.3])
    assert properties.oscillator_strengths == pytest.approx([w0 * 0.4, w1 * 0.8], abs=1e-6)
    if with_ecd:
        assert properties.rotatory_strengths == pytest.approx([w0 * 10, -w1 * 20], abs=1e-4)
    else:
        assert properties.rotatory_strengths is None
    assert properties.vibrational_frequencies == pytest.approx([4400.0, 4400.0])
    assert properties.ir_intensities == pytest.approx([w0 * 100, w1 * 100], abs=1e-3)


def test_boltzmann_weights_require_positive_temperature(monkeypatch: pytest.MonkeyPatch) -> None:
    """Reject a non-positive temperature."""
    result = _ensemble_result(monkeypatch, rowan.ConformerProperties())
    with pytest.raises(ValueError, match="temperature must be positive"):
        result.get_ensemble_properties(0)


def test_broaden_spectrum_conserves_intensity() -> None:
    """Broadened lines integrate to the summed stick intensities."""
    grid = [i * 0.01 for i in range(-1000, 2001)]
    for lineshape in ("gaussian", "lorentzian"):
        spectrum = rowan.broaden_spectrum([3.0, 5.0], [1.0, 0.5], grid, 0.3, lineshape)
        assert sum(spectrum) * 0.01 == pytest.approx(1.5, rel=0.02)

    with pytest.raises(ValueError, match="same length"):
        rowan.broaden_spectrum([3.0], [1.0, 0.5], grid, 0.3)


def test_frequencies_reach_the_optimization_settings(monkeypatch: pytest.MonkeyPatch) -> None:
    """Request frequencies through the final-method preset."""
    client = _mock_submission(monkeypatch)
    rowan.submit_conformer_search_workflow(initial_molecule="CCO", frequencies=True)

    workflow_data = client.post.call_args.kwargs["json"]["workflow_data"]
    assert workflow_data["multistage_opt_settings"]["frequencies"] is True
