"""Test logP workflow submission payloads."""

from unittest.mock import MagicMock

from pytest import MonkeyPatch, mark

import rowan
from rowan.workflows import logp

_PALLADIUM_CHLORIDE = rowan.Molecule.from_xyz("3\n\nPd 0 0 0\nCl 2.3 0 0\nCl -2.3 0 0")


@mark.parametrize(
    ("smiles", "molecule", "input_field", "expected"),
    [
        ("CCO", None, "initial_smiles", "CCO"),
        (
            None,
            _PALLADIUM_CHLORIDE,
            "initial_molecule",
            _PALLADIUM_CHLORIDE.to_stjames().model_dump(mode="json"),
        ),
    ],
)
def test_submit_cosmors_logp_payload(
    monkeypatch: MonkeyPatch,
    smiles: str | None,
    molecule: rowan.Molecule | None,
    input_field: str,
    expected: str | dict[str, object],
) -> None:
    """Send the selected input in both workflow data and the submission envelope."""
    client = MagicMock()
    client.__enter__.return_value = client
    client.post.return_value.json.return_value = {
        "name": "test",
        "uuid": "workflow-uuid",
        "created_at": "2026-09-15T00:00:00Z",
        "object_status": 0,
        "parent_uuid": "folder-uuid",
        "notes": "",
        "starred": False,
        "public": False,
        "object_type": "logp",
        "object_data": {},
        "email_when_complete": False,
        "credits_charged": 0,
        "object_logfile": "",
    }
    monkeypatch.setattr(logp, "api_client", lambda: client)

    rowan.submit_logp_workflow(smiles, method="cosmors", initial_molecule=molecule)

    payload = client.post.call_args.kwargs["json"]
    assert payload[input_field] == payload["workflow_data"][input_field] == expected
    assert {"initial_smiles", "initial_molecule"}.intersection(payload) == {input_field}
