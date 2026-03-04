"""Physical constants and unit conversions for workflows."""

from stjames.data.elements import BRAGG_RADII, ELEMENT_SYMBOL, ISOTOPES, SYMBOL_ELEMENT

# Energy conversions
HARTREE_TO_KCAL = 627.509474
"""Conversion factor from Hartree to kcal/mol."""

HARTREE_TO_EV = 27.211386245988
"""Conversion factor from Hartree to electron volts."""

HARTREE_TO_KJ = 2625.4996394799
"""Conversion factor from Hartree to kJ/mol."""

# Physical constants
BOLTZMANN_HARTREE_PER_K = 3.1668105e-6
"""Boltzmann constant in Hartree/K."""

__all__ = [
    "BOLTZMANN_HARTREE_PER_K",
    "BRAGG_RADII",
    "ELEMENT_SYMBOL",
    "HARTREE_TO_EV",
    "HARTREE_TO_KCAL",
    "HARTREE_TO_KJ",
    "ISOTOPES",
    "SYMBOL_ELEMENT",
]
