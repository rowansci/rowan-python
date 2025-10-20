import asyncio
import copy
import logging
import time
from typing import Iterable, Literal, TypeAlias, TypedDict

import numpy as np
import stjames
from rdkit import Chem
from rdkit.Chem import AllChem

import rowan
from rowan.utils import ATOMIC_NUMBER_TO_ATOMIC_SYMBOL, get_api_key

RdkitMol: TypeAlias = Chem.rdchem.Mol | Chem.rdchem.RWMol
pKaMode = Literal["reckless", "rapid", "careful"]
TautomerMode = Literal["reckless", "rapid", "careful"]
ConformerMode = Literal["reckless", "rapid"]
FAST_METHODS: list[stjames.Method] = [
    *stjames.method.XTB_METHODS,
    *stjames.method.NNP_METHODS,
]


class ConversionError(ValueError):
    pass


class NoConformersError(Exception):
    pass


class MethodTooSlowError(Exception):
    pass


class ChargesResult(TypedDict):
    conformer_index: int
    charges: list[float]


ChargesResults = list[ChargesResult]


class ConformerResult(TypedDict):
    molecule: RdkitMol
    energies: list[float]


class PKaResult(TypedDict):
    element: str
    index: int
    pKa: float


class PKaResults(TypedDict):
    acidic_pkas: list[PKaResult]
    basic_pkas: list[PKaResult]


class TautomerResult(TypedDict):
    molecule: RdkitMol
    predicted_relative_energy: float
    weight: float


TautomerResults = list[TautomerResult]


class ConformerEnergyResult(TypedDict):
    conformer_index: int
    energy: float


class OptimizeResult(TypedDict):
    molecule: RdkitMol
    energies: list[float]


def apply_nest_asyncio() -> None:
    try:
        asyncio.get_running_loop()
    except RuntimeError:
        return
    try:
        import nest_asyncio  # type: ignore [import-untyped]

        nest_asyncio.apply()
    except ImportError:
        pass


# actually apply it
apply_nest_asyncio()


def _get_rdkit_mol_from_uuid(calculation_uuid: str) -> RdkitMol:
    stjames_mol_dict = rowan.retrieve_calculation_molecules(calculation_uuid)[-1]

    return Chem.MolFromXYZBlock(stjames.Molecule(**stjames_mol_dict).to_xyz())


def _embed_rdkit_mol(rdkm: RdkitMol) -> RdkitMol:
    try:
        AllChem.SanitizeMol(rdkm)  # type: ignore [attr-defined]
    except Exception as e:
        raise ValueError("Molecule could not be generated -- invalid chemistry!") from e

    rdkm = AllChem.AddHs(rdkm)  # type: ignore [attr-defined]
    try:
        assert AllChem.EmbedMolecule(rdkm, maxAttempts=200) >= 0  # type: ignore [attr-defined]
    except AssertionError as e:
        status1 = AllChem.EmbedMolecule(rdkm, maxAttempts=200, useRandomCoords=True)  # type: ignore [attr-defined]
        if status1 < 0:
            raise ValueError("Cannot embed molecule!") from e
    try:
        assert AllChem.MMFFOptimizeMolecule(rdkm, maxIters=200) >= 0  # type: ignore [attr-defined]
    except AssertionError:
        pass

    return rdkm


def _rdkit_to_stjames(rdkm: RdkitMol, cid: int = 0) -> stjames.Molecule:
    return stjames.Molecule.from_rdkit(rdkm, cid=cid)


def run_pka(
    mol: RdkitMol,
    mode: pKaMode = "rapid",
    timeout: int = 600,
    name: str = "pKa API Workflow",
    pka_range: tuple[int, int] = (2, 12),
    deprotonate_elements: list[int] | None = None,
    protonate_elements: list[int] | None = None,
    folder_uuid: str | None = None,
) -> PKaResults:
    """
    Calculate the pKa of a Molecule.

    :param mol: RDKit Molecule
    :param mode: pKa calculation Mode. See [list of available modes](https://github.com/rowansci/stjames-public/blob/master/stjames/mode.py)
    for options.
    :param timeout: time in seconds before the Workflow times out
    :param name: name for the job
    :param pka_range: range of pKa values to calculate
    :param deprotonate_elements: elements to deprotonate
    :param protonate_elements: elements to protonate
    :param folder_uuid: folder UUID
    :return: dictionary of pKa values indexed by atom
    """
    protonate_elements = protonate_elements or [7]
    deprotonate_elements = deprotonate_elements or [7, 8, 16]

    return asyncio.run(
        _single_pka(
            mol,
            mode,
            timeout,
            name,
            pka_range,
            deprotonate_elements,
            protonate_elements,
            folder_uuid,
        )
    )


def batch_pka(
    mols: Iterable[RdkitMol],
    mode: pKaMode = "rapid",
    timeout: int = 600,
    name: str = "pKa API Workflow",
    pka_range: tuple[int, int] = (2, 12),
    deprotonate_elements: list[int] | None = None,
    protonate_elements: list[int] | None = None,
    folder_uuid: str | None = None,
) -> list[PKaResults]:
    """
    Calculate the pKa of a batch of Molecules.

    :param mols: list of RDKit Molecules
    :param mode: pKa calculation mode. See [list of available modes](https://github.com/rowansci/stjames-public/blob/master/stjames/mode.py)
    for options.
    :param timeout: time in seconds before the Workflow times out
    :param name: name for the job
    :param pka_range: range of pKa values to calculate
    :param deprotonate_elements: elements to deprotonate
    :param protonate_elements: elements to protonate
    :return: list of dictionary of pKa values indexed by atom
    """
    protonate_elements = protonate_elements or [7]
    deprotonate_elements = deprotonate_elements or [7, 8, 16]

    async def _run():
        tasks = [
            _single_pka(
                mol,
                mode,
                timeout,
                name,
                pka_range,
                deprotonate_elements,
                protonate_elements,
                folder_uuid,
            )
            for mol in mols
        ]
        return await asyncio.gather(*tasks)

    return asyncio.run(_run())


async def _single_pka(
    mol: RdkitMol,
    mode: pKaMode = "rapid",
    timeout: int = 600,
    name: str = "pKa API Workflow",
    pka_range: tuple[int, int] = (2, 12),
    deprotonate_elements: list[int] | None = None,
    protonate_elements: list[int] | None = None,
    folder_uuid: str | None = None,
) -> PKaResults:
    """
    Calculate the pKa of a Molecule.

    :param mol: RDKit Molecule
    :param mode: pKa calculation mode. See [list of available modes](https://github.com/rowansci/stjames-public/blob/master/stjames/mode.py)
    for options.
    :param timeout: time in seconds before the Workflow times out
    :param name: name for the job
    :param pka_range: range of pKa values to calculate
    :param deprotonate_elements: elements to deprotonate
    :param protonate_elements: elements to protonate
    :param folder_uuid: folder UUID
    :return: dictionary of pKa values
    """
    get_api_key()
    protonate_elements = protonate_elements or [7]
    deprotonate_elements = deprotonate_elements or [7, 8, 16]

    post = rowan.submit_workflow(
        name=name,
        workflow_type="pka",
        initial_molecule=_rdkit_to_stjames(mol),
        workflow_data={
            "pka_range": pka_range,
            "deprotonate_elements": deprotonate_elements,
            "deprotonate_atoms": [],
            "protonate_elements": protonate_elements,
            "protonate_atoms": [],
            "mode": mode,
        },
        folder_uuid=folder_uuid,
    )

    start = time.time()
    while not post.is_finished():
        await asyncio.sleep(5)
        if time.time() - start > timeout:
            raise TimeoutError("Workflow timed out")

    data = rowan.retrieve_workflow(post.uuid).data

    if not data:
        raise Exception("Could not retrieve workflow data")

    acidic_pkas: list[PKaResult] = []
    for microstate in data.get("conjugate_bases", []):
        atomic_number = data.get("initial_molecule", {})["atoms"][microstate["atom_index"] - 1][
            "atomic_number"
        ]
        acidic_pkas.append(
            {
                "element": ATOMIC_NUMBER_TO_ATOMIC_SYMBOL[str(atomic_number)],
                "index": microstate["atom_index"],
                "pKa": round(microstate["pka"], 2),
            }
        )

    basic_pkas: list[PKaResult] = []
    for microstate in data.get("conjugate_bases", []):
        atomic_number = data.get("initial_molecule", {})["atoms"][microstate["atom_index"] - 1][
            "atomic_number"
        ]

        basic_pkas.append(
            {
                "element": ATOMIC_NUMBER_TO_ATOMIC_SYMBOL[str(atomic_number)],
                "index": microstate["atom_index"],
                "pKa": round(microstate["pka"], 2),
            }
        )

    return {"acidic_pkas": acidic_pkas, "basic_pkas": basic_pkas}


def run_tautomers(
    mol: RdkitMol,
    mode: TautomerMode = "reckless",
    timeout: int = 600,
    name: str = "Tautomers API Workflow",
    folder_uuid: str | None = None,
) -> TautomerResults:
    """
    Generate possible tautomers of a Molecule.

    :param mol: RDKit Molecule
    :param mode: Tautomer mode. See [list of available modes](https://github.com/rowansci/stjames-public/blob/master/stjames/mode.py)
    for options.
    :param timeout: time in seconds before the Workflow times out
    :param name: name for the job
    :param folder_uuid: folder UUID
    :return: list of dictionaries containing RDKit Molecule, relative energies, and weights
    """
    return asyncio.run(_single_tautomers(mol, mode, timeout, name, folder_uuid))


def batch_tautomers(
    mols: Iterable[RdkitMol],
    mode: TautomerMode = "reckless",
    timeout: int = 600,
    name: str = "Tautomers API Workflow",
    folder_uuid: str | None = None,
) -> list[TautomerResults]:
    """
    Generate possible tautomers of a Molecule.

    :param mols: RDKit Molecule
    :param mode: Tautomer mode. See [list of available modes](https://github.com/rowansci/stjames-public/blob/master/stjames/mode.py)
    for options.
    :param timeout: time in seconds before the Workflow times out
    :param name: name for the job
    :param folder_uuid: folder UUID
    :return: list of lists of dictionaries containing RDKit Molecule, relative energies, and weights
    """

    async def _run():
        tasks = [_single_tautomers(mol, mode, timeout, name, folder_uuid) for mol in mols]
        return await asyncio.gather(*tasks)

    return asyncio.run(_run())


async def _single_tautomers(
    mol: RdkitMol,
    mode: TautomerMode = "reckless",
    timeout: int = 600,
    name: str = "Tautomers API Workflow",
    folder_uuid: str | None = None,
) -> TautomerResults:
    """
    Generate possible tautomers of a Molecule.

    :param mol: RDKit Molecule
    :param mode: Tautomer mode
    :param timeout: time in seconds before the Workflow times out
    :param name: name for the job
    :param folder_uuid: folder UUID
    :return: dictionaries containing RDKit Molecule, relative energy, and weight
    """
    get_api_key()

    post = rowan.submit_workflow(
        name=name,
        workflow_type="tautomers",
        initial_molecule=_rdkit_to_stjames(mol),
        workflow_data={"mode": mode},
        folder_uuid=folder_uuid,
    )

    start = time.time()
    while not post.is_finished():
        await asyncio.sleep(5)
        if time.time() - start > timeout:
            raise TimeoutError("Workflow timed out")

    data = rowan.retrieve_workflow(post.uuid).data

    if not data:
        raise Exception("Could not retrieve workflow data")

    return [
        {
            "molecule": _get_rdkit_mol_from_uuid(tautomer["structures"][0]["uuid"]),
            "predicted_relative_energy": round(tautomer["predicted_relative_energy"], 2),
            "weight": round(tautomer["weight"], 5),
        }
        for tautomer in data.get("tautomers", [])
    ]


def run_energy(
    mol: RdkitMol,
    method: str = "aimnet2_wb97md3",
    engine: str = "aimnet2",
    mode: str = "auto",
    timeout: int = 600,
    name: str = "Energy API Workflow",
    folder_uuid: str | None = None,
) -> list[ConformerEnergyResult]:
    """
    Computes the energy for the given molecule.

    :param mol: RDKit Molecule
    :param method: Method to use for the calculation.
    See [list of available methods](https://github.com/rowansci/stjames-public/blob/master/stjames/method.py)
    :param engine: Engine to run the energy. See [list of available engines](https://github.com/rowansci/stjames-public/blob/master/stjames/engine.py)
    :param mode: Mode to run the energy. See [list of available modes](https://github.com/rowansci/stjames-public/blob/master/stjames/mode.py)
    for options.
    :param timeout: time in seconds before the Workflow times out
    :param name: name for the job
    :param folder_uuid: folder UUID
    :raises: MethodTooSlowError if the method is invalid
    :return: dictionary with the energy in Hartree and the conformer index
    """
    return asyncio.run(_single_energy(mol, method, engine, mode, timeout, name, folder_uuid))


def batch_energy(
    mols: Iterable[RdkitMol],
    method: str = "aimnet2_wb97md3",
    engine: str = "aimnet2",
    mode: str = "auto",
    timeout: int = 600,
    name: str = "Energy API Workflow",
    folder_uuid: str | None = None,
) -> list[list[ConformerEnergyResult]]:
    """
    Computes the energy for the given molecule.

    :param mols: RDKit Molecule
    :param method: Method to use for the calculation.
    See [list of available methods](https://github.com/rowansci/stjames-public/blob/master/stjames/method.py)
    :param engine: Engine to run the energy. See [list of available engines](https://github.com/rowansci/stjames-public/blob/master/stjames/engine.py)
    :param mode: Mode to run the energy. See [list of available modes](https://github.com/rowansci/stjames-public/blob/master/stjames/mode.py)
    for options.
    :param timeout: time in seconds before the Workflow times out
    :param name: name for the job
    :param folder_uuid: folder UUID
    :raises: MethodTooSlowError if the method is invalid
    :return: list of dictionaries with the energy in Hartree and the conformer index
    """

    async def _run():
        tasks = [
            _single_energy(mol, method, engine, mode, timeout, name, folder_uuid) for mol in mols
        ]
        return await asyncio.gather(*tasks)

    return asyncio.run(_run())


async def _single_energy(
    mol: RdkitMol,
    method: str = "aimnet2_wb97md3",
    engine: str = "aimnet2",
    mode: str = "auto",
    timeout: int = 600,
    name: str = "Energy API Workflow",
    folder_uuid: str | None = None,
) -> list[ConformerEnergyResult]:
    """
    Computes the energy for the given molecule.

    :param mol: RDKit Molecule
    :param method: Method to use for the calculation.
    See [list of available methods](https://github.com/rowansci/stjames-public/blob/master/stjames/method.py)
    :param engine: Engine to run the energy. See [list of available engines](https://github.com/rowansci/stjames-public/blob/master/stjames/engine.py)
    :param mode: Mode to run the energy
    :param timeout: time in seconds before the Workflow times out
    :param name: name for the job
    :param folder_uuid: folder UUID
    :raises: MethodTooSlowError if the method is invalid
    :return: dictionary with the energy in Hartree and the conformer index
    """
    get_api_key()
    method = stjames.Method(method)

    if mol.GetNumConformers() == 0:
        mol = _embed_rdkit_mol(mol)
        if mol.GetNumConformers() == 0:
            raise NoConformersError("This molecule has no conformers")

    if method not in FAST_METHODS:
        raise MethodTooSlowError(
            "This method is too slow; try running this through our web interface."
        )

    workflow_uuids = []
    for conformer in mol.GetConformers():
        cid = conformer.GetId()
        stjames_mol = _rdkit_to_stjames(mol, cid)
        post = rowan.submit_workflow(
            name=name,
            workflow_type="basic_calculation",
            initial_molecule=stjames_mol,
            workflow_data={
                "settings": {
                    "method": method.value,
                    "corrections": [],
                    "tasks": ["energy"],
                    "mode": mode,
                    "opt_settings": {"constraints": []},
                },
                "engine": engine,
            },
            folder_uuid=folder_uuid,
        )

        workflow_uuids.append(post.uuid)

    start = time.time()
    while not all(rowan.retrieve_workflow(uuid).is_finished() for uuid in workflow_uuids):
        await asyncio.sleep(5)
        if time.time() - start > timeout:
            raise TimeoutError("Workflow timed out")

    results = [rowan.retrieve_workflow(uuid).data for uuid in workflow_uuids]

    energies = [
        rowan.retrieve_calculation_molecules(data["calculation_uuid"])[-1]["energy"]
        for data in results
        if data is not None
    ]

    return [{"conformer_index": index, "energy": energy} for index, energy in enumerate(energies)]


def run_optimize(
    mol: RdkitMol,
    method: str = "aimnet2_wb97md3",
    engine: str = "aimnet2",
    mode: str = "auto",
    return_energies: bool = False,
    timeout: int = 600,
    name: str = "Optimize API Workflow",
    folder_uuid: str | None = None,
) -> OptimizeResult:
    """
    Optimize each of a molecule's conformers and then return the molecule.

    :param mol: RDKit Molecule
    :param method: Method to use for the calculation.
    See [list of available methods](https://github.com/rowansci/stjames-public/blob/master/stjames/method.py)
    :param engine: Engine to run the optimization. See [list of available engines](https://github.com/rowansci/stjames-public/blob/master/stjames/engine.py)
    :param return_energies: whether to return energies in Hartree too
    :param mode: Mode to run the optimization. See [list of available modes](https://github.com/rowansci/stjames-public/blob/master/stjames/mode.py)
    for options.
    :param timeout: time in seconds before the Workflow times out
    :param name: name for the job
    :param folder_uuid: folder UUID
    :raises: MethodTooSlowError if the method is invalid
    :return: dictionary with the optimized conformer(s) and optional list of energies per conformer
    """
    return asyncio.run(
        _single_optimize(mol, method, engine, mode, return_energies, timeout, name, folder_uuid)
    )


def batch_optimize(
    mols: Iterable[RdkitMol],
    method: str = "aimnet2_wb97md3",
    engine: str = "aimnet2",
    mode: str = "auto",
    return_energies: bool = False,
    timeout: int = 600,
    name: str = "Optimize API Workflow",
    folder_uuid: str | None = None,
) -> list[OptimizeResult]:
    """
    Optimize each of a Molecule's conformers and then return the Molecule.

    :param mols: input Molecule
    :param method: Method to use for the calculation.
    See [list of available methods](https://github.com/rowansci/stjames-public/blob/master/stjames/method.py)
    :param engine: Engine to run the optimization. See [list of available engines](https://github.com/rowansci/stjames-public/blob/master/stjames/engine.py)
    :param mode: Mode to run the optimization. See [list of available modes](https://github.com/rowansci/stjames-public/blob/master/stjames/mode.py)
    for options.
    :param return_energies: whether to return energies in Hartree too
    :param timeout: time in seconds before the Workflow times out
    :param name: name for the job
    :param folder_uuid: folder UUID
    :raises: MethodTooSlowError if the Method is invalid
    :return: dictionaries with optimized conformer(s) and optional list of energies per conformer
    """

    async def _run():
        tasks = [
            _single_optimize(mol, method, engine, mode, return_energies, timeout, name, folder_uuid)
            for mol in mols
        ]
        return await asyncio.gather(*tasks)

    return asyncio.run(_run())


async def _single_optimize(
    mol: RdkitMol,
    method: str = "aimnet2_wb97md3",
    engine: str = "aimnet2",
    mode: str = "auto",
    return_energies: bool = False,
    timeout: int = 600,
    name: str = "Optimize API Workflow",
    folder_uuid: str | None = None,
) -> OptimizeResult:
    """
    Optimize each of a molecule's conformers and then return the molecule.

    :param mol: RDKit Molecule
    :param method: Method to use for the calculation.
    See [list of available methods](https://github.com/rowansci/stjames-public/blob/master/stjames/method.py)
    :param engine: Engine to run the optimization. See [list of available engines](https://github.com/rowansci/stjames-public/blob/master/stjames/engine.py)
    :param mode: Mode to run the optimization. See [list of available modes](https://github.com/rowansci/stjames-public/blob/master/stjames/mode.py)
    :param return_energies: whether to return energies in Hartree too
    :param timeout: time in seconds before the Workflow times out
    :param name: name for the job
    :param folder_uuid: folder UUID
    :raises: MethodTooSlowError if the method is invalid
    :return: dictionary with the optimized conformer(s) and optional list of energies per conformer
    """
    get_api_key()
    method = stjames.Method(method)

    if mol.GetNumConformers() == 0:
        mol = _embed_rdkit_mol(mol)
        if mol.GetNumConformers() == 0:
            raise NoConformersError("This molecule has no conformers")

    if method not in FAST_METHODS:
        raise MethodTooSlowError(
            "This method is too slow; try running this through our web interface."
        )

    optimized_mol = copy.deepcopy(mol)

    workflow_uuids = []
    for conformer in mol.GetConformers():
        cid = conformer.GetId()
        stjames_mol = _rdkit_to_stjames(mol, cid)

        post = rowan.submit_workflow(
            name=name,
            workflow_type="basic_calculation",
            initial_molecule=stjames_mol,
            workflow_data={
                "settings": {
                    "method": method.value,
                    "corrections": [],
                    "tasks": ["optimize"],
                    "mode": mode,
                    "opt_settings": {"constraints": []},
                },
                "engine": engine,
            },
            folder_uuid=folder_uuid,
        )

        workflow_uuids.append(post.uuid)

    start = time.time()
    while not all(rowan.retrieve_workflow(uuid).is_finished() for uuid in workflow_uuids):
        await asyncio.sleep(5)
        if time.time() - start > timeout:
            raise TimeoutError("Workflow timed out")

    results = [rowan.retrieve_workflow(uuid).data for uuid in workflow_uuids]
    calculations = [
        rowan.retrieve_calculation_molecules(data["calculation_uuid"])
        for data in results
        if data is not None
    ]
    optimization_atoms = [cacluation[-1]["atoms"] for cacluation in calculations]
    optimized_positions = [[atom["position"] for atom in atoms] for atoms in optimization_atoms]

    energies = [cacluation[-1]["energy"] for cacluation in calculations]

    for i, conformer in enumerate(optimized_mol.GetConformers()):
        conformer.SetPositions(np.array(optimized_positions[i]))

    return {
        "molecule": mol,
        "energies": energies if return_energies else [],
    }


def run_conformers(
    mol: RdkitMol,
    num_conformers=10,
    method: str = "aimnet2_wb97md3",
    mode: str = "rapid",
    return_energies: bool = False,
    timeout: int = 600,
    name: str = "Conformer API Workflow",
    folder_uuid: str | None = None,
) -> ConformerResult:
    """
    Generate conformers for a Molecule.

    :param mol: RDKit Molecule
    :param num_conformers: number of conformers to generate
    :param method: Method to use for the calculation.
    See [list of available methods](https://github.com/rowansci/stjames-public/blob/master/stjames/method.py)
    :param mode: Mode for conformer generation. See [list of available modes](https://github.com/rowansci/stjames-public/blob/master/stjames/mode.py)
    for options.
    :param return_energies: whether to return energies in Hartree too
    :param timeout: time in seconds before the Workflow times out
    :param name: name for the job
    :param folder_uuid: folder UUID
    :return: dictionary with the RDKit Molecule and energies
    """
    return asyncio.run(
        _single_conformers(
            mol,
            num_conformers,
            method,
            mode,
            return_energies,
            timeout,
            name,
            folder_uuid,
        )
    )


def batch_conformers(
    mols: Iterable[RdkitMol],
    num_conformers=10,
    method: str = "aimnet2_wb97md3",
    mode: str = "rapid",
    return_energies: bool = False,
    timeout: int = 600,
    name: str = "Conformer API Workflow",
    folder_uuid: str | None = None,
) -> list[ConformerResult]:
    """
    Generate conformers for a Molecule.

    :param mols: RDKit molecule object
    :param num_conformers: number of conformers to generate
    :param method: Method to use for the calculation.
    See [list of available methods](https://github.com/rowansci/stjames-public/blob/master/stjames/method.py)
    :param mode: conformer generation mode. See [list of available modes](https://github.com/rowansci/stjames-public/blob/master/stjames/mode.py)
    for options.
    :param return_energies: whether to return energies in Hartree too
    :param timeout: time in seconds before the Workflow times out
    :param name: name for the job
    :return: list of dictonaries with RDKit Molecule and energies
    """

    async def _run():
        tasks = [
            _single_conformers(
                mol,
                num_conformers,
                method,
                mode,
                return_energies,
                timeout,
                name,
                folder_uuid,
            )
            for mol in mols
        ]
        return await asyncio.gather(*tasks)

    return asyncio.run(_run())


async def _single_conformers(
    mol: RdkitMol,
    num_conformers=10,
    method: str = "aimnet2_wb97md3",
    mode: str = "rapid",
    return_energies: bool = False,
    timeout: int = 600,
    name: str = "Conformer API Workflow",
    folder_uuid: str | None = None,
) -> ConformerResult:
    """
    Generate conformers for a molecule.

    :param mol: RDKit Molecule
    :param num_conformers: number of conformers to generate
    :param method: Method to use for the calculation.
    See [list of available methods](https://github.com/rowansci/stjames-public/blob/master/stjames/method.py)
    :param mode: conformer generation mode. See [list of available modes](https://github.com/rowansci/stjames-public/blob/master/stjames/mode.py)
    for options.
    :param return_energies: whether to return energies in Hartree too
    :param timeout: time in seconds before the Workflow times out
    :param name: name for the job
    :param folder_uuid: folder UUID
    :return: dictionary with RDKit molecule and energies
    """
    get_api_key()
    method = stjames.Method(method)

    if mol.GetNumConformers() == 0:
        mol = _embed_rdkit_mol(mol)
        if mol.GetNumConformers() == 0:
            raise NoConformersError("This molecule has no conformers")

    if method not in FAST_METHODS:
        raise MethodTooSlowError(
            "This method is too slow; try running this through our web interface."
        )

    post = rowan.submit_workflow(
        name=name,
        workflow_type="conformer_search",
        initial_molecule=_rdkit_to_stjames(mol),
        workflow_data={
            "conf_gen_mode": "rapid",
            "mode": mode,
            "mso_mode": "manual",
            "multistage_opt_settings": {
                "mode": "manual",
                "optimization_settings": [
                    {
                        "method": method.value,
                        "tasks": ["optimize"],
                        "corrections": [],
                        "mode": "auto",
                    }
                ],
                "solvent": None,
                "transition_state": False,
                "constraints": [],
            },
        },
        folder_uuid=folder_uuid,
    )

    start = time.time()
    while not post.is_finished():
        await asyncio.sleep(5)
        if time.time() - start > timeout:
            raise TimeoutError("Workflow timed out")

    data = rowan.retrieve_workflow(post.uuid).data

    if data is None:
        raise NoConformersError("This molecule has no conformers")

    sorted_data = sorted(
        zip(data["energies"], data["conformer_uuids"], strict=True),
        key=lambda x: x[0],
    )

    if len(sorted_data) < num_conformers:
        logging.warning(
            "Number of conformers requested is greater than number of conformers available"
        )
    num_conformers = min(num_conformers, len(sorted_data))

    # Extract the UUIDs of the lowest n energies
    lowest_n_uuids = [item[1][0] for item in sorted_data[:num_conformers]]
    lowest_energies = [item[0] for item in sorted_data[:num_conformers]]

    AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers)  # type: ignore [attr-defined]

    for i, conformer in enumerate(mol.GetConformers()):
        atoms = rowan.retrieve_calculation_molecules(lowest_n_uuids[i])[-1]["atoms"]
        pos = [atom["position"] for atom in atoms]
        conformer.SetPositions(np.array(pos))

    return {
        "molecule": mol,
        "energies": lowest_energies if return_energies else [],
    }


def run_charges(
    mol: RdkitMol,
    method: str = "aimnet2_wb97md3",
    engine: str = "aimnet2",
    mode: str = "auto",
    timeout: int = 600,
    name: str = "Charges API Workflow",
    folder_uuid: str | None = None,
) -> ChargesResults:
    """
    Computes atom-centered charges for the given Molecule.

    :param mol: RDKit Molecule
    :param method: Method to use for the calculation.
    See [list of available methods](https://github.com/rowansci/stjames-public/blob/master/stjames/method.py)
    :param engine: Engine to run the charges. See [list of available engines](https://github.com/rowansci/stjames-public/blob/master/stjames/engine.py)
    :param mode: The mode to run the calculation in. See [list of available modes](https://github.com/rowansci/stjames-public/blob/master/stjames/mode.py)
    for options.
    :param timeout: timeout in seconds
    :param name: name of the job
    :param folder_uuid: folder UUID
    :raises: MethodTooSlowError if the method is invalid
    :return: dictionary with the charges and the conformer index
    """
    return asyncio.run(_single_charges(mol, method, engine, mode, timeout, name, folder_uuid))


def batch_charges(
    mols: Iterable[RdkitMol],
    method: str = "aimnet2_wb97md3",
    engine: str = "aimnet2",
    mode: str = "auto",
    timeout: int = 600,
    name: str = "Charges API Workflow",
    folder_uuid: str | None = None,
) -> list[ChargesResults]:
    """
    Computes atom-centered charges for the given Molecules.

    :param mols: RDKit Molecule
    :param method: Method to use for the calculation.
    See [list of available methods](https://github.com/rowansci/stjames-public/blob/master/stjames/method.py)
    :param engine: Engine to run the charges. See [list of available engines](https://github.com/rowansci/stjames-public/blob/master/stjames/engine.py)
    :param mode: The mode to run the calculation in. See [list of available modes](https://github.com/rowansci/stjames-public/blob/master/stjames/mode.py)
    for options.
    :param timeout: timeout in seconds
    :param name: name of the job
    :param folder_uuid: folder UUID
    :raises: MethodTooSlowError if the method is invalid
    :return: list of dictionaries with the charges and the conformer index
    """

    async def _run():
        tasks = [
            _single_charges(mol, method, engine, mode, timeout, name, folder_uuid) for mol in mols
        ]
        return await asyncio.gather(*tasks)

    return asyncio.run(_run())


async def _single_charges(
    mol: RdkitMol,
    method: str = "aimnet2_wb97md3",
    engine: str = "aimnet2",
    mode: str = "auto",
    timeout: int = 600,
    name: str = "Energy API Workflow",
    folder_uuid: str | None = None,
) -> ChargesResults:
    """
    Computes atom-centered charges for the given Molecule.

    :param mol: RDKit Molecule
    :param method: Method to use for the calculation.
    See [list of available methods](https://github.com/rowansci/stjames-public/blob/master/stjames/method.py)
    :param engine: Engine to run the charges. See [list of available engines](https://github.com/rowansci/stjames-public/blob/master/stjames/engine.py)
    :param mode: The mode to run the calculation in. See [list of available modes](https://github.com/rowansci/stjames-public/blob/master/stjames/mode.py)
    for options.
    :param timeout: timeout in seconds
    :param name: name of the job
    :param folder_uuid: folder UUID
    :raises: MethodTooSlowError if the method is invalid
    :return: dictionary with the charges and the conformer index
    """
    get_api_key()
    method = stjames.Method(method)

    if mol.GetNumConformers() == 0:
        mol = _embed_rdkit_mol(mol)
        if mol.GetNumConformers() == 0:
            raise NoConformersError("This molecule has no conformers")

    if method not in FAST_METHODS:
        raise MethodTooSlowError(
            "This method is too slow; try running this through our web interface."
        )

    workflow_uuids = []
    for conformer in mol.GetConformers():
        cid = conformer.GetId()

        post = rowan.submit_workflow(
            name=name,
            workflow_type="basic_calculation",
            initial_molecule=_rdkit_to_stjames(mol, cid),
            workflow_data={
                "settings": {
                    "method": method.value,
                    "corrections": [],
                    "tasks": ["charge"],
                    "mode": mode,
                    "opt_settings": {"constraints": []},
                },
                "engine": engine,
            },
            folder_uuid=folder_uuid,
        )

        workflow_uuids.append(post.uuid)

    start = time.time()
    while not all(rowan.retrieve_workflow(uuid).is_finished() for uuid in workflow_uuids):
        await asyncio.sleep(5)
        if time.time() - start > timeout:
            raise TimeoutError("Workflow timed out")

    def grab_charges(uuid: str) -> list[float]:
        """Grab mulliken charges by UUID of workflow."""
        data = rowan.retrieve_workflow(uuid).data
        if data is None:
            raise KeyError("Workflow data not found")
        molecules = rowan.retrieve_calculation_molecules(data["calculation_uuid"])
        return molecules[-1]["mulliken_charges"]

    return [
        {"conformer_index": i, "charges": grab_charges(uuid)}
        for i, uuid in enumerate(workflow_uuids)
    ]
