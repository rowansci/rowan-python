"""Molecule class for representing molecular structures and computed properties."""

from pathlib import Path
from typing import Self

import stjames
from pydantic import BaseModel, ConfigDict, PrivateAttr


class Molecule(BaseModel):
    """
    A molecular structure with optional computed properties.

    Can be created from SMILES, XYZ, or directly from atoms. Wraps stjames.Molecule
    internally but provides a cleaner interface.

    :ivar charge: Molecular charge
    :ivar multiplicity: Spin multiplicity
    :ivar atoms: List of atoms with positions
    :ivar energy: Electronic energy (Hartree)
    :ivar smiles: SMILES string representation
    """

    _stjames: stjames.Molecule = PrivateAttr()

    model_config = ConfigDict(arbitrary_types_allowed=True)

    def __init__(self, _stjames: stjames.Molecule | None = None, **kwargs):
        """Initialize from stjames.Molecule or keyword arguments."""
        super().__init__()
        if _stjames is not None:
            self._stjames = _stjames
        else:
            # Build stjames.Molecule from kwargs
            self._stjames = stjames.Molecule(**kwargs)

    @classmethod
    def from_smiles(cls, smiles: str) -> Self:
        """
        Create a molecule from a SMILES string.

        :param smiles: SMILES string
        :return: Molecule instance
        """
        stj = stjames.Molecule.from_smiles(smiles)
        return cls(_stjames=stj)

    @classmethod
    def from_xyz(cls, xyz_string: str, charge: int = 0, multiplicity: int = 1) -> Self:
        """
        Create a molecule from an XYZ string.

        :param xyz_string: XYZ format string
        :param charge: Molecular charge (default 0)
        :param multiplicity: Spin multiplicity (default 1)
        :return: Molecule instance
        """
        stj = stjames.Molecule.from_xyz(xyz_string)
        # Set charge and multiplicity after creation
        stj = stjames.Molecule(
            atoms=stj.atoms,
            charge=charge,
            multiplicity=multiplicity,
        )
        return cls(_stjames=stj)

    @classmethod
    def from_xyz_file(cls, path: str | Path, charge: int = 0, multiplicity: int = 1) -> Self:
        """
        Create a molecule from an XYZ file.

        :param path: Path to XYZ file
        :param charge: Molecular charge (default 0)
        :param multiplicity: Spin multiplicity (default 1)
        :return: Molecule instance
        """
        xyz_string = Path(path).read_text()
        return cls.from_xyz(xyz_string, charge=charge, multiplicity=multiplicity)

    @classmethod
    def _from_stjames(cls, stj: stjames.Molecule) -> Self:
        """Create from an existing stjames.Molecule (internal use)."""
        return cls(_stjames=stj)

    def __repr__(self) -> str:
        n_atoms = len(self.atoms)
        e_str = f"{self.energy:.6f}" if self.energy is not None else "None"
        return f"<Molecule atoms={n_atoms} energy={e_str}>"

    @property
    def charge(self) -> int:
        """Molecular charge."""
        return self._stjames.charge

    @property
    def multiplicity(self) -> int:
        """Spin multiplicity."""
        return self._stjames.multiplicity

    @property
    def atoms(self) -> list[stjames.Atom]:
        """List of atoms with positions."""
        return self._stjames.atoms

    @property
    def num_atoms(self) -> int:
        """Number of atoms."""
        return len(self._stjames.atoms)

    @property
    def smiles(self) -> str | None:
        """SMILES string representation."""
        return self._stjames.smiles

    @property
    def energy(self) -> float | None:
        """Electronic energy (Hartree)."""
        return self._stjames.energy

    @property
    def charges(self) -> list[float] | None:
        """Partial charges on each atom."""
        return self._stjames.mulliken_charges

    @property
    def spin_densities(self) -> list[float] | None:
        """Spin densities on each atom."""
        return self._stjames.mulliken_spin_densities

    @property
    def dipole(self) -> tuple[float, float, float] | None:
        """Dipole moment vector (Debye)."""
        return self._stjames.dipole

    @property
    def homo_lumo_gap(self) -> float | None:
        """HOMO-LUMO gap (Hartree)."""
        return self._stjames.homo_lumo_gap

    @property
    def vibrational_modes(self) -> list[stjames.VibrationalMode] | None:
        """Vibrational modes from frequency calculation."""
        return self._stjames.vibrational_modes

    @property
    def frequencies(self) -> list[float] | None:
        """Vibrational frequencies (cm^-1)."""
        modes = self._stjames.vibrational_modes
        if modes is None:
            return None
        return [m.frequency for m in modes]

    @property
    def zero_point_energy(self) -> float | None:
        """Zero-point vibrational energy (Hartree)."""
        return self._stjames.zero_point_energy

    @property
    def thermal_energy_correction(self) -> float | None:
        """Thermal correction to energy (Hartree)."""
        return self._stjames.thermal_energy_corr

    @property
    def thermal_enthalpy_correction(self) -> float | None:
        """Thermal correction to enthalpy (Hartree)."""
        return self._stjames.thermal_enthalpy_corr

    @property
    def thermal_free_energy_correction(self) -> float | None:
        """Thermal correction to Gibbs free energy (Hartree)."""
        return self._stjames.thermal_free_energy_corr

    @property
    def gradient(self) -> list[tuple[float, float, float]] | None:
        """Nuclear gradient (Hartree/Bohr)."""
        return self._stjames.gradient

    def to_xyz(self) -> str:
        """Convert to XYZ format string."""
        return self._stjames.to_xyz()

    def to_dict(self) -> dict:
        """Convert to dictionary representation."""
        return self._stjames.model_dump(mode="json")

    def _to_stjames(self) -> stjames.Molecule:
        """Get the underlying stjames.Molecule (internal use)."""
        return self._stjames


__all__ = ["Molecule"]
