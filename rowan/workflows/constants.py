"""Physical constants and unit conversions for workflows."""

import math

# Energy conversions
HARTREE_TO_KCAL = 627.509474

# Physical constants
BOLTZMANN_HARTREE_PER_K = 3.1668105e-6

# Temperature of the ensemble properties Rowan computes, in Kelvin
ENSEMBLE_TEMPERATURE = 298.15


def to_relative_kcal(energies: list[float]) -> list[float]:
    """Convert absolute Hartree energies to relative kcal/mol.

    Args:
        energies: absolute energies in Hartree

    Returns:
        energies relative to the minimum, in kcal/mol. Empty list if input is empty
    """
    if not energies:
        return []
    min_e = min(energies)
    return [(e - min_e) * HARTREE_TO_KCAL for e in energies]


def boltzmann_weights(energies: list[float], temperature: float) -> list[float]:
    """Boltzmann weights of states from their energies.

    Args:
        energies: energies in Hartree
        temperature: temperature in Kelvin

    Returns:
        weight of each state (sum to 1.0)

    Raises:
        ValueError: temperature is not positive or there are no energies

    >>> [round(w, 3) for w in boltzmann_weights([0.0, 0.001], 298.15)]
    [0.743, 0.257]
    """
    if temperature <= 0:
        raise ValueError(f"temperature must be positive, got {temperature} K")
    if not energies:
        raise ValueError("No energies to weight")
    min_e = min(energies)
    kt = BOLTZMANN_HARTREE_PER_K * temperature
    factors = [math.exp(-(e - min_e) / kt) for e in energies]
    total = sum(factors)
    return [f / total for f in factors]
