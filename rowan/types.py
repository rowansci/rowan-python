"""Shared type aliases for the Rowan package."""

from typing import Any, TypeAlias

import stjames
from rdkit import Chem

from .molecule import Molecule as RowanMolecule

RdkitMol: TypeAlias = Chem.rdchem.Mol | Chem.rdchem.RWMol
StJamesMolecule: TypeAlias = stjames.Molecule
MoleculeInput: TypeAlias = dict[str, Any] | RowanMolecule | StJamesMolecule | RdkitMol | str
SolventInput: TypeAlias = stjames.Solvent | str | None
SMILES: TypeAlias = str
UUID: TypeAlias = str
