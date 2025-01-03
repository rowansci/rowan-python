from rdkit import Chem
from rdkit.Chem import AllChem
from typing import Literal, TypeAlias
import rowan
import copy
import math
from rowan.utils import get_api_key
import stjames
import time
import numpy as np
import cctk

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
    def __init__(self, message):
        super().__init__(message)
        self.message = message

class MethodTooSlowError(Exception):
    def __init__(self, message):
        super().__init__(message)
        self.message = message

def _get_rdkit_mol_from_uuid(calculation_uuid: str) -> RdkitMol:
    stjames_mol_dict = rowan.Calculation.retrieve(calculation_uuid)["molecules"][-1]
    stjames_mol = stjames.Molecule(**stjames_mol_dict)
    rdkm = Chem.MolFromXYZBlock(stjames_mol.to_xyz())
    return rdkm

def embed_rdkit_mol(rdkm: RdkitMol):
    try:
        AllChem.SanitizeMol(rdkm)
    except Exception as e:
        raise ValueError(f"Molecule could not be generated -- invalid chemistry!\n{e}")

    rdkm = AllChem.AddHs(rdkm)
    try:
        status1 = AllChem.EmbedMolecule(rdkm, maxAttempts=200)
        assert status1 >= 0
    except Exception as e:
        status1 = AllChem.EmbedMolecule(rdkm, maxAttempts=200, useRandomCoords=True)
        if status1 < 0:
            raise ValueError(f"Cannot embed molecule! Error: {e}")
    try:
        status2 = AllChem.MMFFOptimizeMolecule(rdkm, maxIters=200)
        assert status2 >= 0
    except AssertionError:
        pass
    
    return rdkm

def rdkit_to_cctk(rdkm: RdkitMol, cid: int = 0) -> cctk.Molecule:
    if len(rdkm.GetConformers()) == 0:
        rdkm = embed_rdkit_mol(rdkm)
    try:
        nums = [atom.GetAtomicNum() for atom in rdkm.GetAtoms()]
        geom = rdkm.GetConformers()[cid].GetPositions()
        return cctk.Molecule(nums, geom, charge=Chem.GetFormalCharge(rdkm))
    except IndexError as e:
        raise ConversionError("RDKit molecule does not have a conformer with the given ID") from e
    
def cctk_to_stjames(cmol: cctk.Molecule) -> stjames.Molecule:
    atomic_numbers = cmol.atomic_numbers.view(np.ndarray)
    geometry = cmol.geometry.view(np.ndarray)
    atoms = []
    for i in range(cmol.num_atoms()):
        atoms.append(stjames.Atom(atomic_number=atomic_numbers[i], position=geometry[i]))

    return stjames.Molecule(atoms=atoms, charge=cmol.charge, multiplicity=cmol.multiplicity)

def rdkit_to_stjames(rdkm: RdkitMol, cid: int = 0) -> stjames.Molecule:
    cmol = rdkit_to_cctk(rdkm, cid)
    return cctk_to_stjames(cmol)

# add default ranges / elements
def pka(mol: RdkitMol,
        mode: pKaMode = "rapid",
        timeout: int = 600,
        name: str = "pKa API Workflow",
        pka_range: tuple[int, int] = (2, 12),
        deprotonate_elements: list[int] = [7, 8, 16],
        protonate_elements: list[int] = [7]) -> tuple[dict[int, float], dict[int, float]]:
    """
    Calculate the pKa of a molecule.
    :param mol: RDKit molecule object
    :return: Estimated pKa (float or list, depending on implementation)
    """
    get_api_key()
    post = rowan.Workflow.submit(
        name=name,
        workflow_type="pka",
        initial_molecule=rdkit_to_stjames(mol),
        workflow_data={
            "pka_range": pka_range,
            "deprotonate_elements": deprotonate_elements,
            "deprotonate_atoms": [],
            "protonate_elements": protonate_elements,
            "protonate_atoms": [],
            "mode": mode,
        },
    )
    
    start = time.time()
    while not rowan.Workflow.is_finished(post["uuid"]):
        time.sleep(5)
        if time.time() - start > timeout:
            raise TimeoutError("Workflow timed out")
        
    result = rowan.Workflow.retrieve(post["uuid"])
    print(f"Res: {result}")
    
    acidic_pkas = []
    for microstate in result["object_data"]["conjugate_bases"]:
        symbol = cctk.helper_functions.get_symbol(result["object_data"]["initial_molecule"]["atoms"][microstate["atom_index"]-1]["atomic_number"])
        acidic_pkas.append({
            "element": symbol,
            "index": microstate["atom_index"],
            "pKa": round(microstate["pka"], 2)
        })
        

    basic_pkas = []
    for microstate in result["object_data"]["conjugate_acids"]:
        symbol = cctk.helper_functions.get_symbol(result["object_data"]["initial_molecule"]["atoms"][microstate["atom_index"]-1]["atomic_number"])
        basic_pkas.append({
            "element": symbol,
            "index": microstate["atom_index"],
            "pKa": round(microstate["pka"], 2)
        })

    return {"acidic_pkas": acidic_pkas, "basic_pkas": basic_pkas}


def tautomers(mol: RdkitMol,
              mode: TautomerMode = "reckless",
              timeout: int = 600,
              name: str = "Tautomers API Workflow") -> list[tuple[RdkitMol, float]]:
    """
    Generate possible tautomers of a molecule.
    :param mol: RDKit molecule object
    :return: List of RDKit molecule objects
    """
    get_api_key()
    post = rowan.Workflow.submit(
        name=name,
        workflow_type="tautomers",
        initial_molecule=rdkit_to_stjames(mol),
        workflow_data={
            "mode": mode,
        },
    )
    
    start = time.time()
    while not rowan.Workflow.is_finished(post["uuid"]):
        time.sleep(5)
        if time.time() - start > timeout:
            raise TimeoutError("Workflow timed out")
        
    result = rowan.Workflow.retrieve(post["uuid"])
    
    tautomers = []
    for tautomer in result["object_data"]["tautomers"]:
        rdkit_mol = _get_rdkit_mol_from_uuid(tautomer["structures"][0]["uuid"])
        tautomers.append({
            "molecule": rdkit_mol,
            "predicted_relative_energy": round(tautomer["predicted_relative_energy"], 2),
            "weight": round(tautomer["weight"], 5)
        })
        
    #return relative weights too 
    return tautomers

def energy(
    mol: RdkitMol,
    method: str = "aimnet2_wb97md3",
    timeout: int = 600,
    name: str = "Energy API Workflow",
):
    """
    Computes the energy for the given molecule.

    :param mol: the input molecule
    :param method: the method with which to compute the molecule's energy
    :raises: MethodTooSlowError if the method is invalid
    :returns: the energy, in Hartree
    """
    get_api_key()

    method = stjames.Method(method)

    if mol.GetNumConformers() == 0:
        mol = embed_rdkit_mol(mol)
        if mol.GetNumConformers() == 0:
            raise NoConformersError("This molecule has no conformers")

    if method not in FAST_METHODS:
        raise MethodTooSlowError(
            "This method is too slow; try running this through our web interface."
        )

    workflow_uuids = []
    for conformer in mol.GetConformers():
        cid = conformer.GetId()
        stjames_mol = rdkit_to_stjames(mol, cid)
        get_api_key()
        post = rowan.Workflow.submit(
            name=name,
            workflow_type="basic_calculation",
            initial_molecule=rdkit_to_stjames(mol),
            workflow_data={
            "settings": {
                        "method": "aimnet2_wb97md3",
                        "corrections": [],
                        "tasks": [
                            "energy"
                        ],
                        "mode": "auto",
                        "opt_settings": {
                            "constraints": []
                        }
                    },
            "engine": "aimnet2"
                },
        )
        
        workflow_uuids.append(post["uuid"])
        
    start = time.time()
    while not all(rowan.Workflow.is_finished(uuid) for uuid in workflow_uuids):
        time.sleep(5)
        if time.time() - start > timeout:
            raise TimeoutError("Workflow timed out")
        
    workflow_results = [rowan.Workflow.retrieve(uuid) for uuid in workflow_uuids]
    energies = [rowan.Calculation.retrieve(workflow["object_data"]["calculation_uuid"])["molecules"][-1]["energy"] for workflow in workflow_results]

    return [{"conformer_index": index, "energy": energy} for index, energy in enumerate(energies)]
        
        
def optimize(
    mol: RdkitMol,
    method: str = "aimnet2_wb97md3",
    return_energies: bool = False,
    timeout: int = 600,
    name: str = "Optimize API Workflow",
) -> RdkitMol | tuple[RdkitMol, dict[int, float]]:
    """
    Optimize each of a molecule's conformers and then return the molecule.

    :param mol: the input molecule
    :param method: the method with which to compute the molecule's energy
    :param return_energies: whether to return energies in Hartree too
    :raises: MethodTooSlowError if the method is invalid
    :returns: the molecule, with optimized conformers, and optionally a list of energies per conformer too

    """
    get_api_key()

    method = stjames.Method(method)
    
    if mol.GetNumConformers() == 0:
        mol = embed_rdkit_mol(mol)
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
        stjames_mol = rdkit_to_stjames(mol, cid)
        get_api_key()
        post = rowan.Workflow.submit(
            name=name,
            workflow_type="basic_calculation",
            initial_molecule=rdkit_to_stjames(mol),
            workflow_data={
                "settings": {
                    "method": "aimnet2_wb97md3",
                    "corrections": [],
                    "tasks": [
                        "optimize"
                    ],
                    "mode": "auto",
                    "opt_settings": {
                        "constraints": []
                    }
                },
                "engine": "aimnet2"
            },
        )
        
        workflow_uuids.append(post["uuid"])
        
    start = time.time()
    while not all(rowan.Workflow.is_finished(uuid) for uuid in workflow_uuids):
        time.sleep(5)
        if time.time() - start > timeout:
            raise TimeoutError("Workflow timed out")
        
    workflow_results = [rowan.Workflow.retrieve(uuid) for uuid in workflow_uuids]
    calculations = [rowan.Calculation.retrieve(workflow["object_data"]["calculation_uuid"]) for workflow in workflow_results]
    optimization_atoms = [cacluation["molecules"][-1]["atoms"] for cacluation in calculations]
    optimized_positions = [[atom["position"] for atom in atoms] for atoms in optimization_atoms]
    
    energies = [cacluation["molecules"][-1]["energy"] for cacluation in calculations]
    
    for i, conformer in enumerate(optimized_mol.GetConformers()):
        conformer.SetPositions(np.array(optimized_positions[i]))
        
    return_dict = {"rdkit molecule": mol}
    
    if return_energies:
        return_dict["energies"] = energies

    return return_dict

def conformers(mol: RdkitMol, num_conformers=10,
               method: str = "aimnet2_wb97md3",
                return_energies: bool = False,
                timeout: int = 600,
                name: str = "Conformer API Workflow"):
    """
    Generate conformers for a molecule.
    :param mol: RDKit molecule object
    :param num_conformers: Number of conformers to generate
    :return: List of RDKit molecule objects
    """
    
    get_api_key()

    method = stjames.Method(method)
    
    if mol.GetNumConformers() == 0:
        mol = embed_rdkit_mol(mol)
        if mol.GetNumConformers() == 0:
            raise NoConformersError("This molecule has no conformers")

    if method not in FAST_METHODS:
        raise MethodTooSlowError(
            "This method is too slow; try running this through our web interface."
        )
    
    post = rowan.Workflow.submit(
        name=name,
        workflow_type="conformer_search",
        initial_molecule=rdkit_to_stjames(mol),
        workflow_data={
        "conf_gen_mode": "rapid",
        "mode": "rapid",
        "mso_mode": "manual",
        "multistage_opt_settings": {
            "mode": "manual",
            "optimization_settings": [
                {
                    "method": "aimnet2_wb97md3",
                    "tasks": [
                        "optimize"
                    ],
                    "corrections": [],
                    "mode": "auto"
                }
            ],
            "solvent": None,
            "transition_state": False,
            "constraints": []
        }
    },
    )
    
    start = time.time()
    while not rowan.Workflow.is_finished(post["uuid"]):
        time.sleep(5)
        if time.time() - start > timeout:
            raise TimeoutError("Workflow timed out")
        
    result = rowan.Workflow.retrieve(post["uuid"])
    
    sorted_data = sorted(zip(result["object_data"]["energies"], result["object_data"]["conformer_uuids"]), key=lambda x: x[0])

    # Extract the UUIDs of the lowest n energies
    lowest_n_uuids = [item[1][0] for item in sorted_data[:num_conformers]]
    lowest_energies = [item[0] for item in sorted_data[:num_conformers]]
    print(lowest_energies)
    
    AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers)
    
    for i, conformer in enumerate(mol.GetConformers()):
        atoms = rowan.Calculation.retrieve(lowest_n_uuids[i])["molecules"][-1]["atoms"]
        pos = [atom["position"] for atom in atoms]
        conformer.SetPositions(np.array(pos))
        
    return_dict = {"rdkit molecule": mol}
    
    if return_energies:
        return_dict["energies"] = lowest_energies

    return return_dict
