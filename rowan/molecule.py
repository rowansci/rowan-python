"""Molecule class for representing molecular structures and computed properties."""

import math
import random
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
    def molecules_from_sdf(cls, path: str | Path) -> list[Self]:
        """
        Read all records from an SDF file as a list of molecules.

        Each record becomes one Molecule, in file order, with atom ordering preserved
        across records -- so a multi-conformer SDF reads as a consistent ensemble.

        :param path: path to the SDF file
        :returns: one Molecule per record
        """
        return [cls.from_stjames(mol) for mol in stjames.Molecule.molecules_from_sdf(path)]

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
    def name(self) -> str | None:
        """Optional molecule name."""
        return self._stjames.name

    @name.setter
    def name(self, value: str | None) -> None:
        self._stjames.name = value

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

    @property
    def symmetry(self) -> str | int | None:
        """Space group number (1–230) or point group symbol."""
        return self._stjames.symmetry

    @property
    def xrd_peaks(self) -> list[tuple[int, int, int, float]] | None:
        """XRD reflections as (h, k, l, intensity) tuples."""
        return self._stjames.x_ray_diffraction_peaks

    @property
    def band_structure(self) -> "stjames.BandStructure | None":
        """Electronic band structure and DOS (periodic systems only)."""
        return self._stjames.band_structure

    @property
    def band_gap(self) -> float | None:
        """Band gap (Hartree)."""
        bs = self._stjames.band_structure
        return bs.band_gap if bs is not None else None

    @property
    def density_of_states(self) -> list[tuple[float, float]] | None:
        """Total DOS as (energy in Hartree, k-weighted count) pairs (Fermi level = 0)."""
        bs = self._stjames.band_structure
        return bs.total_density_of_states if bs is not None else None

    @property
    def elastic_tensor(
        self,
    ) -> (
        tuple[
            tuple[float, float, float, float, float, float],
            tuple[float, float, float, float, float, float],
            tuple[float, float, float, float, float, float],
            tuple[float, float, float, float, float, float],
            tuple[float, float, float, float, float, float],
            tuple[float, float, float, float, float, float],
        ]
        | None
    ):
        """Elastic stiffness matrix in GPa (Voigt order: xx yy zz yz xz xy)."""
        return self._stjames.elastic_tensor

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

    def perturb(self, stddev: float = 0.005) -> "Molecule":
        """
        Return a copy with random Gaussian noise added to each atom position.

        Useful for breaking symmetry before resubmitting an optimization (e.g. when
        a calculation is stuck in a saddle point or converges to an unwanted geometry).

        :param stddev: standard deviation of the Gaussian displacement (Angstroms)
        :returns: new Molecule with perturbed coordinates
        """

        def _gauss() -> float:
            u1 = random.random()
            u2 = random.random()
            return math.sqrt(-2.0 * math.log(u1)) * math.cos(2.0 * math.pi * u2) * stddev

        new_atoms = [
            stjames.Atom(
                atomic_number=atom.atomic_number,
                position=(
                    atom.position[0] + _gauss(),
                    atom.position[1] + _gauss(),
                    atom.position[2] + _gauss(),
                ),
                mass=atom.mass,
            )
            for atom in self._stjames.atoms
        ]
        new_mol = self._stjames.model_copy(update={"atoms": new_atoms})
        return Molecule(_stjames=new_mol)

    def displace_along_mode(self, mode: stjames.VibrationalMode, displacement: float) -> "Molecule":
        """
        Return a copy with atom positions displaced along a vibrational normal mode.

        Requires a prior frequency calculation. For a transition state, the imaginary
        mode has a negative frequency and points toward the reactant or product.

        :param mode: vibrational mode to displace along, from ``vibrational_modes``
        :param displacement: displacement distance along the normalized mode (Angstroms)
        :returns: new Molecule with displaced coordinates
        :raises ValueError: if the mode has no displacements or zero-norm displacements
        """
        raw = mode.displacements
        if not raw:
            raise ValueError(
                "Vibrational mode has no displacements. "
                "Rerun the calculation with frequencies=True to displace along a mode."
            )
        norm = math.sqrt(sum(x**2 + y**2 + z**2 for x, y, z in raw))
        if norm == 0:
            raise ValueError(
                "Vibrational mode has zero-norm displacements. "
                "This indicates malformed frequency output, try displacing along a different mode."
            )

        new_atoms = [
            stjames.Atom(
                atomic_number=atom.atomic_number,
                position=(
                    atom.position[0] + (dx / norm) * displacement,
                    atom.position[1] + (dy / norm) * displacement,
                    atom.position[2] + (dz / norm) * displacement,
                ),
                mass=atom.mass,
            )
            for atom, (dx, dy, dz) in zip(self._stjames.atoms, raw, strict=True)
        ]
        new_mol = self._stjames.model_copy(update={"atoms": new_atoms})
        return Molecule(_stjames=new_mol)


def load_named_ligands(path: Path | str) -> dict[str, Molecule]:
    """
    Load named ligands from a multi-molecule file.

    Molecule names are read from the title field of each record. Use this when
    ligand identity needs to be preserved - for example, building a ligand dict
    for an RBFE workflow where names are used as keys throughout submission and
    results. Every record must therefore have a unique name; duplicates raise
    rather than silently collapsing into one ligand.

    Supported formats (all carry per-molecule name fields):
    - SDF / MOL (``.sdf``, ``.mol``) - name from the title line
    - MOL2 (``.mol2``) - name from the ``@<TRIPOS>MOLECULE`` block

    :param path: Path to an SDF, MOL, or MOL2 file.
    :returns: Dict mapping ligand name to Molecule, in file order.
    :raises ValueError: If no valid molecules are found, the format is
        unsupported, or two records share a name.
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

    # require every record to have a unique name rather than dropping data
    names = [name for name, _ in pairs]
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise ValueError(
            f"Ligand names must be unique, but {path} repeats: {', '.join(duplicates)}. "
            "Rename the duplicate records before loading."
        )

    return {name: Molecule(_stjames=stjames.Molecule.from_rdkit(rdkm)) for name, rdkm in pairs}
