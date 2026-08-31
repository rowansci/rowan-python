"""Regression tests for cofolding workflow updates."""

from collections.abc import Iterator
from contextlib import contextmanager
from unittest.mock import MagicMock

import pytest
import stjames
from pydantic import ValidationError
from pytest import MonkeyPatch

import rowan
from rowan.workflows import protein_cofolding as cofolding_module


def _workflow_response() -> dict[str, object]:
    """Build a minimal workflow API response."""
    return {
        "name": "test",
        "uuid": "workflow-uuid",
        "created_at": "2026-08-24T00:00:00Z",
        "object_status": 0,
        "parent_uuid": "folder-uuid",
        "notes": "",
        "starred": False,
        "public": False,
        "object_type": "protein_cofolding",
        "object_data": {},
        "email_when_complete": False,
        "credits_charged": 0,
        "object_logfile": "",
    }


def _mock_submission(monkeypatch: MonkeyPatch) -> MagicMock:
    """Replace the cofolding API client and return its mock client."""
    client = MagicMock()
    client.post.return_value.json.return_value = _workflow_response()

    @contextmanager
    def mock_api_client() -> Iterator[MagicMock]:
        yield client

    monkeypatch.setattr(cofolding_module, "api_client", mock_api_client)
    return client


def test_cofolding_types_are_available_from_rowan() -> None:
    """Keep user-facing sequence, modification, and constraint types under rowan."""
    assert rowan.ProteinSequence is stjames.ProteinSequence
    assert rowan.DNASequence is stjames.DNASequence
    assert rowan.RNASequence is stjames.RNASequence
    assert rowan.ResidueModification is stjames.ResidueModification
    assert rowan.NucleotideModification is stjames.NucleotideModification
    assert rowan.BondConstraint is stjames.BondConstraint


def test_cofolding_typed_sequences_and_bond_constraints(monkeypatch: MonkeyPatch) -> None:
    """Serialize modifications and a covalent constraint through the public API."""
    client = _mock_submission(monkeypatch)
    protein = rowan.ProteinSequence(
        sequence="ASA",
        modifications=[rowan.ResidueModification(position=1, ccd="SEP")],
    )
    dna = rowan.DNASequence(
        sequence="AT",
        modifications=[rowan.NucleotideModification(position=0, ccd="5MC")],
    )
    rna = rowan.RNASequence(
        sequence="AU",
        modifications=[rowan.NucleotideModification(position=1, ccd="PSU")],
    )
    protein_atom = rowan.ConstraintTarget(
        input_type="protein", input_index=0, token_index=1, atom_name="OG"
    )
    ligand_atom = rowan.ConstraintTarget(input_type="ligand", input_index=0, token_index=0)
    bond = rowan.BondConstraint(atom_1=protein_atom, atom_2=ligand_atom)
    rowan.submit_protein_cofolding_workflow(
        initial_protein_sequences=[protein],
        initial_dna_sequences=[dna],
        initial_rna_sequences=[rna],
        initial_smiles_list=["C"],
        bond_constraints=[bond],
    )
    payload = client.post.call_args.kwargs["json"]["workflow_data"]
    assert payload["initial_protein_sequences"][0]["modifications"][0]["ccd"] == "SEP"
    assert payload["initial_dna_sequences"][0]["modifications"][0]["ccd"] == "5MC"
    assert payload["initial_rna_sequences"][0]["modifications"][0]["ccd"] == "PSU"
    assert payload["bond_constraints"][0]["atom_1"]["atom_name"] == "OG"


def test_cofolding_decaf_model_is_serialized(monkeypatch: MonkeyPatch) -> None:
    """Serialize the DeCAF model through the public API."""
    client = _mock_submission(monkeypatch)

    rowan.submit_protein_cofolding_workflow(
        initial_protein_sequences=["ACD"], model=rowan.CofoldingModel.DECAF_BOLTZ
    )

    assert client.post.call_args.kwargs["json"]["workflow_data"]["model"] == "decaf_boltz"


def test_cofolding_bond_constraints_reject_pose_refinement(
    monkeypatch: MonkeyPatch,
) -> None:
    """Delegate incompatible bond-constraint settings to StJames validation."""
    _mock_submission(monkeypatch)
    bond = rowan.BondConstraint(
        atom_1=rowan.ConstraintTarget(
            input_type="protein", input_index=0, token_index=1, atom_name="OG"
        ),
        atom_2=rowan.ConstraintTarget(input_type="ligand", input_index=0, token_index=0),
    )

    with pytest.raises(ValidationError, match="bond_constraints"):
        rowan.submit_protein_cofolding_workflow(
            initial_protein_sequences=["ACD"],
            initial_smiles_list=["C"],
            bond_constraints=[bond],
            do_pose_refinement=True,
        )
