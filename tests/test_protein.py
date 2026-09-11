"""Tests for protein helpers."""

from pathlib import Path
from unittest.mock import Mock, patch

import pytest
from stjames.pdb import pdb_from_mmcif_filestring, pdb_from_pdb_filestring

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


def test_download_mmcif_preserves_large_identifiers(tmp_path: Path) -> None:
    """Load listed protein metadata before exporting overflow-sized identifiers."""
    structure = pdb_from_pdb_filestring(
        "HETATM    1  O   HOH W   1       1.000   2.000   3.000  1.00  0.00           O\nEND\n"
    )
    water = structure.models[0].water.pop("W.1")
    water.atoms[100_000] = water.atoms.pop(1)
    structure.models[0].water["W.10000"] = water
    data = structure.model_dump(mode="json")
    protein = Protein(uuid="protein-uuid", data={"title": "Metadata from list_proteins"})

    def refresh(workflow_uuid: str | None = None) -> None:
        assert workflow_uuid == "workflow-uuid"
        protein.data = data

    with patch.object(Protein, "refresh", side_effect=refresh):
        protein.download_structure(path=tmp_path / "structures", workflow_uuid="workflow-uuid")
    exported = pdb_from_mmcif_filestring((tmp_path / "structures" / "protein-uuid.cif").read_text())
    record = exported.models[0]._atom_records()[0]
    assert record.serial == 100_000
    assert record.residue_number == "10000"
    assert (record.atom.x, record.atom.y, record.atom.z) == (1.0, 2.0, 3.0)
    assert protein.data == data
