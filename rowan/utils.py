import os
import cctk
import stjames
import numpy as np
from contextlib import contextmanager
import httpx

import rowan

from .constants import API_URL


def get_api_key() -> str:
    api_key = os.environ.get("ROWAN_API_KEY")
    if api_key is not None:
        return api_key
    elif hasattr(rowan, "api_key"):
        return rowan.api_key
    else:
        raise ValueError(
            "No API key provided. You can set your API key using 'rowan.api_key = <API-KEY>', or you can set the environment variable ROWAN_API_KEY=<API-KEY>)."
        )


def cctk_to_stjames(molecule: cctk.Molecule) -> stjames.Molecule:
    atomic_numbers = molecule.atomic_numbers.view(np.ndarray)
    geometry = molecule.geometry.view(np.ndarray)

    atoms = list()
    for i in range(molecule.num_atoms()):
        atoms.append(
            stjames.Atom(atomic_number=atomic_numbers[i], position=geometry[i])
        )

    return stjames.Molecule(
        atoms=atoms,
        charge=molecule.charge,
        multiplicity=molecule.multiplicity,
    )


def smiles_to_stjames(smiles: str) -> stjames.Molecule:
    cmol = cctk.Molecule.new_from_smiles(smiles)
    return cctk_to_stjames(cmol)


@contextmanager
def api_client():
    """Wraps `httpx.Client` with Rowan-specific kwargs."""
    with httpx.Client(
        base_url=API_URL,
        headers={"X-API-Key": get_api_key()},
    ) as client:
        yield client
