"""
Collection of classes to hold atomic data and functions to extract and manipulate it
used in BAC related scripts
"""
from __future__ import print_function

from bac_amber_utils import *

class AtomInfo():
    """
    Maps charges to atoms using names, PDB index and RDKit index

    Attributes:
        bound (list): Bound atom indices
        rings (list): Index of rings for which the atom is a member
                      in the list of rings created by RDKit
        name (str): Atom name
        rdkit_name (str): Atom name with surrounding spaces as recorded in RDKit
        charge (float): Charge on atom
    """

    def __init__(self):

        self.bound = []
        self.rings = []
        self.name = ''
        self.rdkit_name = ''
        self.charge = 0.0

        return


class PdbRdkitChargeMap():
    """
    Maps charges to atoms using names, PDB index and RDKit index

    Attributes:
        atom_charges (dict): keys = atom names (str), values = charges (float)
        pdb_idx_name_map (dict): keys = PDB atom indices (int),
                                 values = atom names (str)
        rdkit_idx_name_map (dict): keys = RDKit, 0 based, atom indices (int),
                                   values = atom names (str)
    """

    def __init__(self):

        self.atom_charges = {}
        self.pdb_idx_name_map = {}
        self.rdkit_idx_name_map = {}

        return


    def charge_pdb_idx(self, idx):
        """
        Get charge for selected atom

        Args:
            idx: PDB atom index

        Returns:
            float: atom charge
        """

        name = self.pdb_idx_name_map[idx]

        return self.atom_charges[name]


    def charge_rdkit_idx(self, idx):
        """
        Get charge for selected atom

        Args:
            idx: RDKit atom index

        Returns:
            float: atom charge
        """

        name = self.rdkit_idx_name_map[idx]

        return self.atom_charges[name]


    def name_pdb_idx(self,idx):
        """
        Get name for selected atom

        Args:
            idx: PDB atom index

        Returns:
            str: atom name
        """

        return self.pdb_idx_name_map[idx]


    def name_rdkit_idx(self, idx):
        """
        Get name for selected atom

        Args:
            idx: RDKit atom index

        Returns:
            str: atom name
        """

        return self.rdkit_idx_name_map[idx]


    def create_atom_charge_map(self, names, charges):
        """
        Create a mapping between the input atom names and corresponding charges

        Updates self.atom_charges

        Args:
            names (list): Atom names
            charges (list): Atomic charges

        Returns:

        """

        if len(names) == len(charges):
            self.atom_charges = dict(zip(names,charges))
        else:
            err = "Passed lists are of different lengths"
            raise ValueError(err)

        return

def get_bond_list(mol):
    """
    Get list of the bonds (described by start and end atom indices) in mol

    Args:
        mol (rdkit.Chem.rdchem.Mol):  Information on atomic contents of molecule

    Returns:
        list:  List of pairs of atom indices involved in each bond
    """

    bonds = mol.GetBonds()

    bond_list = []

    for bond in bonds:

        atm1 = bond.GetBeginAtom().GetIdx()
        atm2 = bond.GetEndAtom().GetIdx()

        bond_list.append([atm1,atm2])

    return bond_list


def get_atom_info_rdkit(mol, charge_map=None):
    """
    Extract atom information from RDKit molecule

    Associates atom index (RDKit style) to name, bond, ring and optionally
    charge information.

    Args:
        mol (rdkit.Chem.rdchem.Mol):  Information on atomic contents of molecule

    Kwargs:
        charge_dict (PdbRdkitChargeMap):  Map of atom names and indices to charges

    Returns:
        dict:  keys = atom indices, values = AtomInfo describing atom
    """

    atom_info = {}

    #  Get list of lists of atom indices involved in each ring in the molecule
    rings = mol.GetRingInfo().AtomRings()

    for idx in range(mol.GetNumAtoms()):

        atom_info[idx] = AtomInfo()

        info = atom_info[idx]

        atom = mol.GetAtomWithIdx(idx)

        # RDKit names contain surrounding spaces
        info.rdkit_name = atom.GetPDBResidueInfo().GetName()

        info.name = info.rdkit_name.strip()

        if charge_map:

            info.charge = charge_map.atom_charges[info.name]

        for ring_no in range(len(rings)):

            if idx in rings[ring_no]:

                info.rings.append(ring_no)

    #  Create list of atoms bonded to each atom
    bond_list = get_bond_list(mol)

    for idx1,idx2 in bond_list:

        atom_info[idx1].bound.append(idx2)
        atom_info[idx2].bound.append(idx1)

    return atom_info

def create_charge_idx_map(ac_filename, structure):

    charge_map = PdbRdkitChargeMap()

    charge_map.atom_charges = get_resp_charge_per_atom_from_ac(ac_filename, structure)

    for idx in range(structure.natoms()):

        name = structure.name()[idx]
        pdb_idx = structure.index()[idx]

        charge_map.pdb_idx_name_map[pdb_idx] = name
        charge_map.rdkit_idx_name_map[idx] = name

    return charge_map