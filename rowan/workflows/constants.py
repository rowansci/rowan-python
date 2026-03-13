"""Physical constants and unit conversions for workflows."""

# Energy conversions
HARTREE_TO_KCAL = 627.509474

# Physical constants
BOLTZMANN_HARTREE_PER_K = 3.1668105e-6


def to_relative_kcal(energies: list[float]) -> list[float]:
    """Convert absolute Hartree energies to relative kcal/mol.

    :param energies: Absolute energies in Hartree.
    :returns: Energies relative to the minimum, in kcal/mol. Empty list if input is empty.
    """
    if not energies:
        return []
    min_e = min(energies)
    return [(e - min_e) * HARTREE_TO_KCAL for e in energies]
