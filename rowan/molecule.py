"""Molecule class for representing molecular structures and computed properties."""

from pathlib import Path
from typing import Any, Self

import stjames
from pydantic import BaseModel, ConfigDict, PrivateAttr
from rdkit import Chem


class Molecule(BaseModel):
    """
    Molecular structure with optional computed properties.

    Can be created from SMILES, XYZ, or directly from atoms. Wraps stjames.Molecule
    internally but provides a cleaner interface.

    :param charge: Molecular charge.
    :param multiplicity: Spin multiplicity.
    :param atoms: List of atoms with positions.
    :param energy: Electronic energy (Hartree).
    :param smiles: SMILES string representation.
    """

    _stjames: stjames.Molecule = PrivateAttr()

    model_config = ConfigDict(arbitrary_types_allowed=True)

    def __init__(self, _stjames: stjames.Molecule | None = None, **kwargs):
        """Initialize from stjames.Molecule or keyword arguments."""
        super().__init__()
        if _stjames is not None:
            self._stjames = _stjames
        else:
            self._stjames = stjames.Molecule(**kwargs)

    def __repr__(self) -> str:
        n_atoms = len(self.atoms)
        e_str = f"{self.energy:.6f}" if self.energy is not None else "None"
        return f"<Molecule atoms={n_atoms} energy={e_str}>"

    # -- Constructors (from_*) and converters (to_*) --

    @classmethod
    def from_smiles(cls, smiles: str) -> Self:
        """
        Create molecule from SMILES string.

        :param smiles: SMILES string.
        :returns: Molecule instance.
        """
        return cls(_stjames=stjames.Molecule.from_smiles(smiles))

    @classmethod
    def from_xyz(
        cls, xyz_string: str, charge: int | None = None, multiplicity: int | None = None
    ) -> Self:
        """
        Create Molecule from XYZ string.

        :param xyz_string: XYZ format string
        :param charge: charge
        :param multiplicity: spin multiplicity
        :returns: Molecule
        """
        return cls(
            _stjames=stjames.Molecule.from_xyz(
                xyz_string,
                charge=charge,
                multiplicity=multiplicity,
            )
        )

    @classmethod
    def from_xyz_file(
        cls, path: str | Path, charge: int | None = None, multiplicity: int | None = None
    ) -> Self:
        """
        Create molecule from XYZ file.

        :param path: path to XYZ file
        :param charge: charge
        :param multiplicity: spin multiplicity
        :returns: Molecule
        """
        return cls.from_xyz(Path(path).read_text(), charge=charge, multiplicity=multiplicity)

    @classmethod
    def from_stjames(cls, stj: stjames.Molecule) -> Self:
        """
        Create from stjames.Molecule.

        :param stj: stjames.Molecule
        :returns: Molecule
        """
        return cls(_stjames=stj)

    @classmethod
    def from_atoms(
        cls,
        atoms: list[stjames.Atom],
        charge: int = 0,
        multiplicity: int = 1,
        cell: stjames.PeriodicCell | None = None,
    ) -> Self:
        """
        Create molecule from a list of atoms.

        :param atoms: List of Atom objects.
        :param charge: Molecular charge (default 0).
        :param multiplicity: Spin multiplicity (default 1).
        :param cell: PeriodicCell for periodic boundary conditions.
        :returns: Molecule instance.
        """
        return cls(
            _stjames=stjames.Molecule(
                atoms=atoms,
                charge=charge,
                multiplicity=multiplicity,
                cell=cell,
            )
        )

    def to_xyz(self) -> str:
        """Convert to XYZ format string."""
        return self._stjames.to_xyz()

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary representation."""
        return self._stjames.model_dump(mode="json")

    def to_stjames(self) -> stjames.Molecule:
        """Convert to stjames.Molecule."""
        return self._stjames

    # -- Properties --

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
    def cell(self) -> stjames.PeriodicCell | None:
        """Unit cell for periodic boundary conditions."""
        return self._stjames.cell

    @property
    def gradient(self) -> list[tuple[float, float, float]] | None:
        """Nuclear gradient (Hartree/Bohr)."""
        return self._stjames.gradient

    # -- Geometric utilities --

    def distance(self, i: int, j: int) -> float:
        """
        Calculate distance between two atoms.

        :param i: First atom index (1-based).
        :param j: Second atom index (1-based).
        :returns: Distance in Angstroms.
        """
        return self._stjames.distance(i, j)

    def angle(self, i: int, j: int, k: int, degrees: bool = True) -> float:
        """
        Calculate angle between three atoms (i-j-k).

        :param i: First atom index (1-based).
        :param j: Central atom index (1-based).
        :param k: Third atom index (1-based).
        :param degrees: Return angle in degrees (default) or radians.
        :returns: Angle in degrees or radians.
        """
        return self._stjames.angle(i, j, k, degrees=degrees)

    def dihedral(self, i: int, j: int, k: int, l: int, degrees: bool = True) -> float:
        """
        Calculate dihedral angle between four atoms (i-j-k-l).

        :param i: First atom index (1-based).
        :param j: Second atom index (1-based).
        :param k: Third atom index (1-based).
        :param l: Fourth atom index (1-based).
        :param degrees: Return angle in degrees (default) or radians.
        :returns: Dihedral angle in degrees (0-360) or radians.
        """
        return self._stjames.dihedral(i, j, k, l, degrees=degrees)


def load_named_ligands(path: Path | str) -> dict[str, Molecule]:
    """
    Load named ligands from a multi-molecule file.

    Molecule names are read from the title field of each record. Use this when
    ligand identity needs to be preserved - for example, building a ligand dict
    for an RBFE workflow where names are used as keys throughout submission and
    results.

    Supported formats (all carry per-molecule name fields):
    - SDF / MOL (``.sdf``, ``.mol``) - name from the title line
    - MOL2 (``.mol2``) - name from the ``@<TRIPOS>MOLECULE`` block

    :param path: Path to an SDF, MOL, or MOL2 file.
    :returns: Dict mapping ligand name to Molecule, in file order.
    :raises ValueError: If no valid molecules are found or the format is unsupported.
    """
    path = Path(path)
    suffix = path.suffix.lower()

    if suffix in {".sdf", ".mol"}:
        supplier = Chem.SDMolSupplier(str(path), removeHs=False)
        pairs = [(rdkm.GetProp("_Name"), rdkm) for rdkm in supplier if rdkm is not None]
    elif suffix == ".mol2":
        text = path.read_text()
        blocks = [b for b in text.split("@<TRIPOS>MOLECULE") if b.strip()]
        pairs = []
        for block in blocks:
            rdkm = Chem.MolFromMol2Block("@<TRIPOS>MOLECULE" + block, removeHs=False)
            if rdkm is not None:
                name = rdkm.GetProp("_Name") if rdkm.HasProp("_Name") else ""
                pairs.append((name, rdkm))
    else:
        raise ValueError(f"Unsupported file format: {suffix!r} (expected .sdf, .mol, or .mol2)")

    if not pairs:
        raise ValueError(f"No valid molecules found in {path}")

    return {name: Molecule(_stjames=stjames.Molecule.from_rdkit(rdkm)) for name, rdkm in pairs}
