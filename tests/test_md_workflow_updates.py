"""Regression tests for molecular dynamics workflow updates."""

import struct
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from unittest.mock import MagicMock

import stjames
from pytest import MonkeyPatch
from stjames.pdb import pdb_from_pdb_filestring

import rowan
from rowan.protein import Protein
from rowan.workflows import pose_analysis_md as pose_analysis_md_module
from rowan.workflows import protein_md as protein_md_module
from rowan.workflows.protein_md import ProteinMDResult


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
        "object_type": "protein_md",
        "object_data": {},
        "email_when_complete": False,
        "credits_charged": 0,
        "object_logfile": "",
    }


def _mock_submission(monkeypatch: MonkeyPatch, module: object) -> MagicMock:
    """Replace a workflow module's API client and return its mock client."""
    client = MagicMock()
    client.post.return_value.json.return_value = _workflow_response()

    @contextmanager
    def mock_api_client() -> Iterator[MagicMock]:
        yield client

    monkeypatch.setattr(module, "api_client", mock_api_client)
    return client


def test_md_forcefield_types_are_available_from_rowan() -> None:
    """Keep user-facing force-field types under rowan."""
    assert rowan.ProteinForceField is stjames.ProteinForceField
    assert rowan.WaterForceField is stjames.WaterForceField


def test_md_trajectory_constructors_preserve_their_previous_signatures() -> None:
    """Default newly added fields when constructing trajectory results directly."""
    protein_trajectory = protein_md_module.ProteinMDTrajectory(
        uuid="trajectory-uuid",
        sasa=[],
        polar_sasa=[],
        isotropic_radius_of_gyration=[],
        cluster_centroid_indices=[],
        cluster_indices_by_frame=[],
    )
    pose_trajectory = pose_analysis_md_module.TrajectoryResult(
        uuid="trajectory-uuid",
        ligand_rmsd=[],
        contacts=[],
        sasa=[],
        polar_sasa=[],
        isotropic_radius_of_gyration=[],
        cluster_centroid_indices=[],
        cluster_indices_by_frame=[],
    )

    assert protein_trajectory.binder_rmsd == []
    assert protein_trajectory.mmgbsa_scores == []
    assert protein_trajectory.protein_rmsd == []
    assert protein_trajectory.rmsf == []
    assert protein_trajectory.potential_energy == []
    assert protein_trajectory.mean_structure_uuid is None
    assert protein_trajectory.median_structure_frame_index is None
    assert pose_trajectory.protein_rmsd == []
    assert pose_trajectory.rmsf == []
    assert pose_trajectory.potential_energy == []
    assert pose_trajectory.mmgbsa_scores == []
    assert pose_trajectory.mean_structure_uuid is None
    assert pose_trajectory.median_structure_frame_index is None


def test_protein_md_uses_new_defaults_and_forcefields(monkeypatch: MonkeyPatch) -> None:
    """Use the updated MD defaults and serialize public force-field types."""
    client = _mock_submission(monkeypatch, protein_md_module)
    protein_md_module.submit_protein_md_workflow(
        "protein-uuid",
        protein_ff=rowan.ProteinForceField.FF19SB,
        water_ff=rowan.WaterForceField.OPC,
        validate_forcefield=False,
    )
    payload = client.post.call_args.kwargs["json"]["workflow_data"]
    assert payload["protein_ff"] == "ff19sb"
    assert payload["water_ff"] == "opc"
    assert payload["equilibration_time_ns"] == 0.5
    assert payload["timestep_fs"] == 4
    assert payload["hydrogen_mass"] == 3
    assert payload["water_buffer"] == 8


def test_protein_md_exposes_binder_and_representative_structure_results(
    monkeypatch: MonkeyPatch,
) -> None:
    """Expose binder analyses and cache fetched mean structures."""

    workflow = stjames.ProteinMolecularDynamicsWorkflow(
        protein="protein-uuid",
        trajectories=[
            {
                "uuid": "trajectory-uuid",
                "protein_rmsd": [0.0, 0.4],
                "rmsf": [0.2, 0.3],
                "potential_energy": [-1.1, -1.2],
                "binder_rmsd": [1.0, 1.2],
                "mmgbsa_scores": [-20.0, None],
                "mean_structure_uuid": "mean-uuid",
                "median_structure_frame_index": 7,
            }
        ],
    )
    result = ProteinMDResult(
        workflow_data=workflow.model_dump(mode="json"),
        workflow_type="protein_md",
        workflow_uuid="workflow-uuid",
    )
    trajectory = result.trajectories[0]
    assert trajectory.protein_rmsd == [0.0, 0.4]
    assert trajectory.rmsf == [0.2, 0.3]
    assert trajectory.potential_energy == [-1.1, -1.2]
    assert trajectory.binder_rmsd == [1.0, 1.2]
    assert trajectory.mmgbsa_scores == [-20.0, None]
    assert trajectory.mean_structure_uuid == "mean-uuid"
    assert trajectory.median_structure_frame_index == 7

    mean_protein = Protein(uuid="mean-uuid", name="mean")
    retrieve_mean = MagicMock(return_value=mean_protein)
    monkeypatch.setattr("rowan.workflows._molecular_dynamics.retrieve_protein", retrieve_mean)
    assert result.get_mean_structure() is mean_protein
    assert result.get_mean_structure() is mean_protein
    retrieve_mean.assert_called_once_with("mean-uuid", workflow_uuid="workflow-uuid")


def test_protein_md_accepts_results_without_new_trajectory_fields() -> None:
    """Keep older Protein MD results readable when new fields are absent."""

    old_workflow = stjames.ProteinMolecularDynamicsWorkflow(
        protein="protein-uuid", trajectories=[{"uuid": "old-trajectory"}]
    )
    old_result = ProteinMDResult(
        workflow_data=old_workflow.model_dump(mode="json", exclude_defaults=True),
        workflow_type="protein_md",
        workflow_uuid="old-workflow",
    )
    assert old_result.trajectories[0].binder_rmsd == []
    assert old_result.trajectories[0].mean_structure_uuid is None


def test_protein_md_accepts_legacy_binder_schema() -> None:
    """Keep workflows created with binder-contained SMILES readable."""

    legacy_result = ProteinMDResult(
        workflow_data={
            "protein": "protein-uuid",
            "binder": {"small_molecules": {"LIG": "CC"}},
            "trajectories": [{"uuid": "legacy-trajectory"}],
        },
        workflow_type="protein_md",
        workflow_uuid="legacy-workflow",
    )
    assert legacy_result._workflow.binder.small_molecule_residues == ["LIG"]
    assert legacy_result._workflow.small_molecules == {"LIG": "CC"}


def test_protein_md_serializes_general_binders(monkeypatch: MonkeyPatch) -> None:
    """Serialize protein, small-molecule, and combined binders."""
    client = _mock_submission(monkeypatch, protein_md_module)

    binders: list[tuple[stjames.Binder, dict[str | int, str | None] | None]] = [
        (stjames.Binder(chain_ids=["B", "C"]), None),
        (
            stjames.Binder(small_molecule_residues=["LIG", "INH"]),
            {"LIG": "CC", "INH": "CCC"},
        ),
        (
            stjames.Binder(chain_ids=["B"], small_molecule_residues=["LIG"]),
            {"LIG": "CC"},
        ),
    ]
    for binder, small_molecules in binders:
        protein_md_module.submit_protein_md_workflow(
            "protein-uuid",
            binder=binder,
            small_molecules=small_molecules,
            validate_forcefield=False,
        )
        workflow_payload = client.post.call_args.kwargs["json"]["workflow_data"]
        assert workflow_payload["binder"] == binder.model_dump(mode="json")
        assert workflow_payload.get("small_molecules") == small_molecules


def test_pose_analysis_md_uses_new_defaults_and_forcefields(monkeypatch: MonkeyPatch) -> None:
    """Apply the updated controls and defaults to pose-analysis MD."""
    client = _mock_submission(monkeypatch, pose_analysis_md_module)
    pose_analysis_md_module.submit_pose_analysis_md_workflow(
        "protein-uuid",
        "CC",
        protein_ff=rowan.ProteinForceField.FF19SB,
        water_ff=rowan.WaterForceField.OPC,
        validate_forcefield=False,
    )
    payload = client.post.call_args.kwargs["json"]["workflow_data"]
    assert payload["protein_ff"] == "ff19sb"
    assert payload["water_ff"] == "opc"
    assert payload["equilibration_time_ns"] == 0.5
    assert payload["timestep_fs"] == 4
    assert payload["hydrogen_mass"] == 3
    assert payload["water_buffer"] == 8


def test_pose_analysis_md_exposes_trajectory_analysis_results() -> None:
    """Expose trajectory analyses and representative structures for new and old results."""

    workflow = stjames.PoseAnalysisMolecularDynamicsWorkflow(
        protein="protein-uuid",
        initial_smiles="CC",
        trajectories=[
            {
                "uuid": "trajectory-uuid",
                "protein_rmsd": [0.0, 0.5],
                "rmsf": [0.1, 0.2],
                "potential_energy": [-1.3, -1.4],
                "mmgbsa_scores": [-20.0, None],
                "mean_structure_uuid": "mean-uuid",
                "median_structure_frame_index": 3,
            }
        ],
    )
    result = pose_analysis_md_module.PoseAnalysisMDResult(
        workflow_data=workflow.model_dump(mode="json"),
        workflow_type="pose_analysis_md",
        workflow_uuid="workflow-uuid",
    )
    trajectory = result.trajectories[0]
    assert trajectory.protein_rmsd == [0.0, 0.5]
    assert trajectory.rmsf == [0.1, 0.2]
    assert trajectory.potential_energy == [-1.3, -1.4]
    assert trajectory.mmgbsa_scores == [-20.0, None]
    assert trajectory.mean_structure_uuid == "mean-uuid"
    assert trajectory.median_structure_frame_index == 3

    old_workflow = stjames.PoseAnalysisMolecularDynamicsWorkflow(
        protein="protein-uuid",
        initial_smiles="CC",
        trajectories=[{"uuid": "old-trajectory"}],
    )
    old_result = pose_analysis_md_module.PoseAnalysisMDResult(
        workflow_data=old_workflow.model_dump(mode="json", exclude_defaults=True),
        workflow_type="pose_analysis_md",
        workflow_uuid="old-workflow",
    )
    assert old_result.trajectories[0].protein_rmsd == []
    assert old_result.trajectories[0].rmsf == []
    assert old_result.trajectories[0].potential_energy == []
    assert old_result.trajectories[0].mmgbsa_scores == []
    assert old_result.trajectories[0].mean_structure_uuid is None
    assert old_result.trajectories[0].median_structure_frame_index is None


def test_medoid_structure_download_combines_frame_and_topology(
    monkeypatch: MonkeyPatch, tmp_path: Path
) -> None:
    """Write one compressed-stream frame as PDB without an MDAnalysis dependency."""
    pdb = pdb_from_pdb_filestring(
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N  \n"
        "ATOM      2  CA  ALA A   1       1.000   1.000   1.000  1.00  0.00           C  \n"
        "TER\nEND\n"
    )
    workflow = stjames.ProteinMolecularDynamicsWorkflow(
        protein="protein-uuid",
        minimized_protein_uuid="minimized-uuid",
        trajectories=[{"uuid": "trajectory-uuid", "median_structure_frame_index": 4}],
    )
    result = ProteinMDResult(
        workflow_data=workflow.model_dump(mode="json"),
        workflow_type="protein_md",
        workflow_uuid="workflow-uuid",
    )
    result._cache["minimized_protein"] = Protein(
        uuid="minimized-uuid", name="minimized", data=pdb.model_dump(mode="json")
    )

    response = MagicMock()
    response.headers = {"X-Box-Size-Bytes": "48", "X-Num-Atoms": "2"}
    response.content = bytes(48) + struct.pack("<6d", 2.0, 3.0, 4.0, 5.0, 6.0, 7.0)
    client = MagicMock()
    client.get.return_value = response

    @contextmanager
    def mock_api_client() -> Iterator[MagicMock]:
        yield client

    monkeypatch.setattr("rowan.workflows._molecular_dynamics.api_client", mock_api_client)
    path = result.download_medoid_structure(0, tmp_path, "medoid")
    assert path is not None
    atom_lines = [line for line in path.read_text().splitlines() if line.startswith("ATOM")]
    assert atom_lines[0][30:54] == "   2.000   3.000   4.000"
    assert atom_lines[1][30:54] == "   5.000   6.000   7.000"
    client.get.assert_called_once_with(
        "/trajectory/workflow-uuid/compressed_stream",
        params={"replicate": 0, "start_frame": 4, "num_frames": 1},
    )
