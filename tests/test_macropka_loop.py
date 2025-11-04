"""Submit 5k macropka jobs, grouped into a Rowan folder."""

import os
import time
from pathlib import Path

from dotenv import load_dotenv

PROFILE_COUNT = 5000
PROGRESS_INTERVAL = 500
TARGET_FOLDER_NAME = "batch macropka test"
MARAVIROC_SMILES = (
    "CC(C)NC[C@@H](COc1ccc(cc1Cl)C2=NC(=O)N(C(=N2)C)C3CC3)N4CCN(CC4)C5=NC6=CC=CC=C6N=C5"
)


def load_env() -> None:
    script_path = Path(__file__).resolve()
    candidates = [script_path.parent / ".env", script_path.parent.parent / ".env"]

    for candidate in candidates:
        if load_dotenv(dotenv_path=candidate, override=False):
            break

    if "API_URL" in os.environ and "ROWAN_API_URL" not in os.environ:
        os.environ["ROWAN_API_URL"] = os.environ["API_URL"]


def ensure_folder(rowan_module) -> str:
    for folder in rowan_module.list_folders(size=200):
        if folder.name == TARGET_FOLDER_NAME:
            print(f"Using existing Rowan folder '{folder.name}' ({folder.uuid})")
            return folder.uuid

    folder = rowan_module.create_folder(name=TARGET_FOLDER_NAME)
    print(f"Created Rowan folder '{folder.name}' ({folder.uuid})")
    return folder.uuid


def main() -> None:
    load_env()

    import rowan

    if not (getattr(rowan, "api_key", None) or os.environ.get("ROWAN_API_KEY")):
        raise RuntimeError("Set rowan.api_key or ROWAN_API_KEY before running this script.")

    folder_uuid = ensure_folder(rowan)
    start = time.perf_counter()

    for idx in range(1, PROFILE_COUNT + 1):
        rowan.submit_macropka_workflow(
            initial_smiles=MARAVIROC_SMILES,
            name=f"macropka-maraviroc-{idx}",
            compute_aqueous_solubility=False,
            compute_solvation_energy=False,
            folder_uuid=folder_uuid,
        )

        if idx % PROGRESS_INTERVAL == 0:
            elapsed = time.perf_counter() - start
            print(f"Submitted {idx}/{PROFILE_COUNT} workflows in {elapsed:.1f}s")

    total = time.perf_counter() - start
    print(f"Done. {PROFILE_COUNT} submissions -> folder {folder_uuid}. Total {total:.1f}s")


if __name__ == "__main__":
    main()
