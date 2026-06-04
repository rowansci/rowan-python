"""Shared type aliases for the Rowan package."""

from typing import TypeAlias

import stjames
from rdkit import Chem

from .molecule import Molecule as RowanMolecule

RdkitMol: TypeAlias = Chem.rdchem.Mol | Chem.rdchem.RWMol
StJamesMolecule: TypeAlias = stjames.Molecule
# 3D-structure input: molecule objects carrying coordinates. Excludes bare SMILES strings
# (no geometry) and dicts (internal serialization only).
StructureInput: TypeAlias = RowanMolecule | StJamesMolecule | RdkitMol
SolventInput: TypeAlias = stjames.Solvent | str | None
SMILES: TypeAlias = str
UUID: TypeAlias = str
