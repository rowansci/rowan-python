"""Tests for protein helpers."""

from unittest.mock import Mock, patch

import pytest

from rowan.protein import Protein


def test_get_atom_index_delegates_to_pdb() -> None:
    """Delegate atom lookup to the stjames PDB model."""
    data: dict[str, object] = {"models": []}
    pdb = Mock()
    pdb.get_atom_index.return_value = 42

    with patch("rowan.protein.PDB.model_validate", return_value=pdb) as model_validate:
        index = Protein(uuid="protein-uuid", data=data).get_atom_index(
            "A",
            "701A",
            "C1",
            entity_type="non_polymer",
            model_index=1,
        )

    assert index == 42
    model_validate.assert_called_once_with(data)
    pdb.get_atom_index.assert_called_once_with(
        "A",
        "701A",
        "C1",
        entity_type="non_polymer",
        model_index=1,
    )


def test_get_atom_index_requires_loaded_data() -> None:
    """Reject atom lookup when protein data has not been loaded."""
    with pytest.raises(ValueError, match="Protein data not loaded"):
        Protein(uuid="protein-uuid").get_atom_index("A", 101, "SG")
