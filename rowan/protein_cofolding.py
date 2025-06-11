from typing import Any

import stjames

from .utils import api_client


def submit_protein_cofolding_workflow(
        initial_protein_sequences: list[str],
        initial_smiles_list: list[str] | None = None,
        ligand_binding_affinity_index: int | None = None,
        use_msa_server: bool = True,
        use_potenitals: bool = False,
        name: str = "Cofolding Workflow",
        model: str = stjames.CofoldingModel.BOLTZ_2.value,
        folder_uuid: stjames.UUID | None = None,
    ) -> dict[str, Any]:
        workflow_data = {
            "use_msa_server": use_msa_server,
            "use_potentials": use_potenitals,
            "model": model,
            "ligand_binding_affinity_index": ligand_binding_affinity_index,
            "initial_smiles_list": initial_smiles_list,
            "initial_protein_sequences": initial_protein_sequences
        }
        data = {
            "name": name,
            "folder_uuid": folder_uuid,
            "workflow_type": "protein_cofolding",
            "workflow_data": workflow_data,
        }

        with api_client() as client:
            response = client.post("/workflow", json=data)
            response.raise_for_status()
            return response.json()
