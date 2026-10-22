"""Broaden stick spectra into continuous line shapes."""

import math
from collections.abc import Sequence
from typing import Literal

HARTREE_TO_EV = 27.211386245988

_GAUSSIAN_FWHM_TO_SIGMA = 1 / (2 * math.sqrt(2 * math.log(2)))


def broaden_spectrum(
    positions: Sequence[float],
    intensities: Sequence[float],
    grid: Sequence[float],
    fwhm: float,
    lineshape: Literal["gaussian", "lorentzian"] = "gaussian",
) -> list[float]:
    """Sum line shapes centered on each stick, each with area equal to its intensity.

    Use it on the stick spectra of `ConformerSearchResult.ensemble_properties`, for
    example excitation energies with oscillator or rotatory strengths for UV-vis or ECD,
    or vibrational frequencies with IR intensities. `grid` and `fwhm` share the units of
    `positions`; convert excitation energies from Hartree with `HARTREE_TO_EV` first to
    broaden in eV.

    Args:
        positions: line positions
        intensities: line intensities, parallel to `positions`
        grid: points at which to evaluate the spectrum
        fwhm: full width at half maximum of each line
        lineshape: line shape of each stick

    Returns:
        spectrum at each point of `grid`

    Raises:
        ValueError: `fwhm` is not positive or `intensities` differs in length from `positions`

    >>> [round(y, 4) for y in broaden_spectrum([0.0], [1.0], [0.0, 0.5], fwhm=1.0)]
    [0.9394, 0.4697]
    >>> [round(y, 4) for y in broaden_spectrum([0.0], [2.0], [0.0], 1.0, "lorentzian")]
    [1.2732]
    """
    if fwhm <= 0:
        raise ValueError(f"fwhm must be positive, got {fwhm}")
    if len(positions) != len(intensities):
        raise ValueError(
            f"positions ({len(positions)}) and intensities ({len(intensities)}) "
            "must have the same length"
        )

    if lineshape == "gaussian":
        sigma = fwhm * _GAUSSIAN_FWHM_TO_SIGMA
        norm = 1 / (sigma * math.sqrt(2 * math.pi))

        def shape(dx: float) -> float:
            return norm * math.exp(-0.5 * (dx / sigma) ** 2)

    else:
        gamma = fwhm / 2

        def shape(dx: float) -> float:
            return gamma / (math.pi * (dx**2 + gamma**2))

    return [
        sum(i * shape(x - p) for p, i in zip(positions, intensities, strict=True)) for x in grid
    ]
