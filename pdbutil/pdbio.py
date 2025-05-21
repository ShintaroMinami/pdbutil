import sys
from io import StringIO
from typing import Dict
from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser


MAX_FILE_PATH = 255
DEFAULT_BB_ATOMS = ['N', 'CA', 'C', 'O']
BB_ATOMS_TO_INDEX = {a: i for i, a in enumerate(DEFAULT_BB_ATOMS)}
THREE2ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}
ONE2THREE = {v: k for k, v in THREE2ONE.items()}

_parser= PDBParser(QUIET=True)


def read_pdb(
        pdb: str,
        model_num: int=0,
        bb_atoms_dict: Dict={a: i for i, a in enumerate(DEFAULT_BB_ATOMS)}  # Default mapping,
        ) -> dict:
    """
    Read a PDB file/string and return a dictionary with the following keys:
    
    Inputs:
    pdb: str
        Path to the PDB file or string content of the PDB file.
    model_num: int
        Model number to extract from the PDB file. Default is 0.
    bb_atoms_dict: dict
        Dictionary mapping backbone atom names to their indices.
        Default is {'N':0, 'CA':1, 'C':2, 'O':3}.
    
    Returns:
    dict:
        A dictionary containing the following keys:
    - 'xyz_ca': numpy array of CA atom coordinates
    - 'xyz_bb': numpy array of backbone coordinates
    - 'xyz_aa': list of numpy array of all atom coordinates
    - 'chain': numpy array of chain identifiers
    - 'resnum': numpy array of residue numbers
    - 'res1': numpy array of one letter residue types
    - 'res3': numpy array of three letter residue types
    - 'bfactor': numpy array of B-factors
    - 'insertion': numpy array of insertion codes
    - 'pdbstring': string of the PDB file content
    """
    
    pdb_string = Path(pdb).read_text() if Path(pdb[:MAX_FILE_PATH]).is_file() else pdb
    
    structure = _parser.get_structure('protein', StringIO(pdb_string))

    if len(structure) == 0:
        raise ValueError(f"Error: Invalid PDB file or pdb_string.")
    
    model = structure[model_num]
    
    xyz_ca, xyz_backbone, xyz_allatom = [], [], []
    chains, resnumber, residue1, residue3, bfactor, insertion = [], [], [], [], [], []

    for chain in model:
        chain_id = chain.get_id()
        for residue in chain:
            res3 = residue.get_resname()
            res1 = THREE2ONE.get(res3, 'X')
            _, resnum, ins = residue.get_id()
            residue1.append(res1)
            residue3.append(res3)
            chains.append(chain_id)
            resnumber.append(int(resnum))
            insertion.append(ins)
            xyz_bb = [None] * len(bb_atoms_dict)
            xyz_aa = []
            for atom in residue:
                xyz = atom.get_coord()
                xyz_aa.append(xyz)
                if atom.name == 'CA':
                    xyz_ca.append(xyz)
                    bfac = atom.get_bfactor()
                if atom.name in bb_atoms_dict:
                    xyz_bb[bb_atoms_dict[atom.name]] = xyz
            if any([xyz is None for xyz in xyz_bb]):
                sys.stderr.write(f"Warning: Missing backbone atom in residue {res3} {chain_id}{resnum} and skipped\n")
                continue
            xyz_backbone.append(np.stack(xyz_bb))
            xyz_allatom.append(np.stack(xyz_aa))
            bfactor.append(bfac)

    return {
        'xyz_ca': np.stack(xyz_ca),
        'xyz_bb': np.stack(xyz_backbone),
        'xyz_aa': xyz_allatom,
        'chain': np.array(chains),
        'resnum': np.array(resnumber),
        'res1': np.array(residue1),
        'res3': np.array(residue3),
        'bfactor': np.array(bfactor),
        'insertion': np.array(insertion),
        'pdbstring': pdb_string,
    }


def renumbering_for_each_chain(chain: np.ndarray) -> np.ndarray:
    resnum = np.zeros(len(chain), dtype=int)
    for c in np.unique(chain):
        mask = (chain == c).astype(int)
        cumsum = np.cumsum(mask) * mask
        resnum += cumsum
    return resnum


def write_pdb(
        xyz: np.ndarray=None,
        xyz_bb: np.ndarray=None,
        res3: np.ndarray=None,
        resnum: np.ndarray=None,
        chain: np.ndarray=None,
        bfactor: np.ndarray=None,
        default_resname: str='XXX',
        default_chain: str='A',
        default_bfactor: float=0.0,
        default_insertion: str=' ',
        renumber: bool=False,
        pdb_file: str=None,
        **args,
        ) -> str:
    """
    Write a PDB format string from the given parameters.
    """
    assert (xyz is not None) or (xyz_bb is not None), "xyz or xyz_bb must be provided"
    if xyz is None:
        xyz = xyz_bb
    if type(xyz) != np.ndarray:
        raise TypeError(f"xyz should be a numpy array, got {type(xyz)}")
    if xyz.ndim != 3:
        raise ValueError(f"xyz should be a 3D array, got {xyz.ndim}D")
    atoms = [DEFAULT_BB_ATOMS] * xyz.shape[0]
    res3 = np.full((xyz.shape[0],), default_resname) if res3 is None else res3
    chain = np.full((xyz.shape[0],), default_chain) if chain is None else chain
    bfactor = np.full((xyz.shape[0],), default_bfactor) if bfactor is None else bfactor
    resnum = renumbering_for_each_chain(chain) if (resnum is None) or renumber else resnum
    if xyz.shape[0] != res3.shape[0]:
        raise ValueError(f"xyz and res3 must have the same length, got {xyz.shape[0]} and {res3.shape[0]}")
    if xyz.shape[0] != resnum.shape[0]:
        raise ValueError(f"xyz and resnum must have the same length, got {xyz.shape[0]} and {resnum.shape[0]}")
    if xyz.shape[0] != chain.shape[0]:
        raise ValueError(f"xyz and chain must have the same length, got {xyz.shape[0]} and {chain.shape[0]}")
    if xyz.shape[0] != bfactor.shape[0]:
        raise ValueError(f"xyz and bfactor must have the same length, got {xyz.shape[0]} and {bfactor.shape[0]}")

    insertion = np.full((xyz.shape[0],), default_insertion) # insertion code is not used

    pdb_string, iatom, chain_old = "", 0, None
    for xbb, atom4, r3, i, c, ins, bfac in zip(xyz, atoms, res3, resnum, chain, insertion, bfactor):
        if (chain_old is not None) and (c != chain_old):
            pdb_string += f"TER\n"
        for x, a in zip(xbb, atom4):
            iatom += 1
            pdb_string += f"ATOM  {iatom:5d}  {a:<2}  {r3} {c}{i:4d}{ins:<1}   {x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}{bfac:6.2f}\n"
        chain_old = c
    pdb_string += f"END\n"

    # Write PDB file if a file path is provided
    if pdb_file is not None:
        with open(pdb_file, 'w') as f:
            f.write(pdb_string)
    
    # Return PDB string
    return pdb_string


