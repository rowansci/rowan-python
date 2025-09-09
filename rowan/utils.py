import os
from contextlib import contextmanager
from typing import Generator

import httpx
import stjames

import rowan

from .constants import API_URL


def get_api_key() -> str:
    """
    Get the API key from the environment variable ROWAN_API_KEY or the module-level attribute
    rowan.api_key.

    If neither of these are set, raise a ValueError with a helpful message.

    :return: The API key.
    """
    if hasattr(rowan, "api_key") and rowan.api_key:
        return rowan.api_key
    elif (api_key := os.environ.get("ROWAN_API_KEY")) is not None:
        return api_key

    raise ValueError(
        "No API key provided. You can set your API key using 'rowan.api_key = <API-KEY>',"
        + " or you can set the environment variable ROWAN_API_KEY=<API-KEY>)."
    )


def smiles_to_stjames(smiles: str) -> stjames.Molecule:
    """
    Convert a SMILES string to a `stjames.Molecule` object.

    :param smiles: A string representing the SMILES notation of the molecule.
    :return: A `stjames.Molecule` object created from the given SMILES string.
    """

    return stjames.Molecule.from_smiles(smiles)


@contextmanager
def api_client() -> Generator[httpx.Client, None, None]:
    """Wraps `httpx.Client` with Rowan-specific kwargs."""
    with httpx.Client(
        base_url=API_URL,
        headers={"X-API-Key": get_api_key()},
        timeout=30,
    ) as client:
        yield client


ATOMIC_NUMBER_TO_ATOMIC_SYMBOL = {
    "1": "H",
    "2": "He",
    "3": "Li",
    "4": "Be",
    "5": "B",
    "6": "C",
    "7": "N",
    "8": "O",
    "9": "F",
    "10": "Ne",
    "11": "Na",
    "12": "Mg",
    "13": "Al",
    "14": "Si",
    "15": "P",
    "16": "S",
    "17": "Cl",
    "18": "Ar",
    "19": "K",
    "20": "Ca",
    "21": "Sc",
    "22": "Ti",
    "23": "V",
    "24": "Cr",
    "25": "Mn",
    "26": "Fe",
    "27": "Co",
    "28": "Ni",
    "29": "Cu",
    "30": "Zn",
    "31": "Ga",
    "32": "Ge",
    "33": "As",
    "34": "Se",
    "35": "Br",
    "36": "Kr",
    "37": "Rb",
    "38": "Sr",
    "39": "Y",
    "40": "Zr",
    "41": "Nb",
    "42": "Mo",
    "43": "Tc",
    "44": "Ru",
    "45": "Rh",
    "46": "Pd",
    "47": "Ag",
    "48": "Cd",
    "49": "In",
    "50": "Sn",
    "51": "Sb",
    "52": "Te",
    "53": "I",
    "54": "Xe",
    "55": "Cs",
    "56": "Ba",
    "57": "La",
    "58": "Ce",
    "59": "Pr",
    "60": "Nd",
    "61": "Pm",
    "62": "Sm",
    "63": "Eu",
    "64": "Gd",
    "65": "Tb",
    "66": "Dy",
    "67": "Ho",
    "68": "Er",
    "69": "Tm",
    "70": "Yb",
    "71": "Lu",
    "72": "Hf",
    "73": "Ta",
    "74": "W",
    "75": "Re",
    "76": "Os",
    "77": "Ir",
    "78": "Pt",
    "79": "Au",
    "80": "Hg",
    "81": "Tl",
    "82": "Pb",
    "83": "Bi",
    "84": "Po",
    "85": "At",
    "86": "Rn",
    "87": "Fr",
    "88": "Ra",
    "89": "Ac",
    "90": "Th",
    "91": "Pa",
    "92": "U",
    "93": "Np",
    "94": "Pu",
    "95": "Am",
    "96": "Cm",
    "97": "Bk",
    "98": "Cf",
    "99": "Es",
    "100": "Fm",
    "101": "Md",
    "102": "No",
    "103": "Lr",
    "104": "Rf",
    "105": "Db",
    "106": "Sg",
    "107": "Bh",
    "108": "Hs",
    "109": "Mt",
    "110": "Ds",
    "111": "Rg",
    "112": "Cn",
    "113": "Nh",
    "114": "Fl",
    "115": "Mc",
    "116": "Lv",
    "117": "Ts",
    "118": "Og",
}
