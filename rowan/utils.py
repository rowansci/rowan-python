import os
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass
from typing import Generator

import httpx
import stjames

import rowan

from .constants import API_URL


@dataclass(frozen=True, slots=True)
class _APIContext:
    """Store credentials isolated to the current execution context."""

    api_key: str
    project_uuid: str | None


_api_context: ContextVar[_APIContext | None] = ContextVar("rowan_api_context", default=None)


@contextmanager
def api_credentials(api_key: str, project_uuid: str | None = None) -> Generator[None, None, None]:
    """Temporarily use Rowan credentials in the current execution context.

    Context-local credentials take precedence over module-level and environment configuration.
    Nested contexts restore the previous credentials when they exit, and concurrent threads or
    asynchronous tasks remain isolated from one another.

    :param api_key: Rowan API key
    :param project_uuid: active project UUID, if any
    :yields: control while the credentials are active
    :raises ValueError: API key is empty
    """
    if not api_key:
        raise ValueError("API key cannot be empty.")

    token = _api_context.set(_APIContext(api_key=api_key, project_uuid=project_uuid))
    try:
        yield
    finally:
        _api_context.reset(token)


def get_api_key() -> str:
    """
    Get the API key from the environment variable ROWAN_API_KEY or the module-level attribute
    rowan.api_key.

    If neither of these are set, raise a ValueError with a helpful message.

    :returns: API key.
    """
    if (context := _api_context.get()) is not None:
        return context.api_key
    if hasattr(rowan, "api_key") and rowan.api_key:
        return rowan.api_key
    elif (api_key := os.environ.get("ROWAN_API_KEY")) is not None:
        return api_key

    raise ValueError(
        "No API key provided. You can set your API key using 'rowan.api_key = <API-KEY>',"
        + " or you can set the environment variable ROWAN_API_KEY=<API-KEY>)."
    )


def get_project_uuid() -> str | None:
    """
    Get the active project UUID from the module-level attribute rowan.project_uuid.

    :returns: Project UUID string, or None if not set.
    """
    if (context := _api_context.get()) is not None:
        return context.project_uuid
    if hasattr(rowan, "project_uuid") and rowan.project_uuid:
        return rowan.project_uuid
    return None


def smiles_to_stjames(smiles: str) -> stjames.Molecule:
    """
    Convert a SMILES string to a `stjames.Molecule` object.

    :param smiles: String representing the SMILES notation of the molecule.
    :returns: `stjames.Molecule` object created from the given SMILES string.
    """

    return stjames.Molecule.from_smiles(smiles)


def _raise_for_status(response: httpx.Response) -> None:
    """Response hook that raises HTTPStatusError with the API's detail message if available."""
    try:
        response.raise_for_status()
    except httpx.HTTPStatusError as e:
        try:
            response.read()
            detail = response.json().get("detail")
        except Exception:
            detail = None
        if detail:
            raise httpx.HTTPStatusError(
                f"{e.response.status_code} {detail}",
                request=e.request,
                response=e.response,
            ) from None
        raise


@contextmanager
def api_client() -> Generator[httpx.Client, None, None]:
    """Wraps `httpx.Client` with Rowan-specific kwargs."""
    with httpx.Client(
        base_url=API_URL,
        headers={"X-API-Key": get_api_key()},
        timeout=120,
        event_hooks={"response": [_raise_for_status]},
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
