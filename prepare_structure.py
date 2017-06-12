import shutil
import os
import glob
import subprocess
import sys
import copy
from collections import Counter
from collections import OrderedDict

from rdkit import Chem

from user_selection import *
from bac_ties.bac_amber_utils import read_pdb_prep_pair


class MolInfo:
    """
    Collects molecular information and associated file names.

     Attributes:
        mol (rdkit.Chem.rdchem.Mol):
        struct (sasmol,SasMol):
        atom_types (dict):
        output_basename (str):
        ac_filename (str):
        prep_filename (str):
        frcmod_filename (str):
        pdb_filename (str):
        original_pdb (str):
        map_names_to_old_pdb (dict):

    """

    def __init__(self):
        """
        Provides empty attributes.
        """

        self.mol = None
        self.struct = None
        self.atom_types = {}
        self.output_basename = ''
        self.ac_filename = ''
        self.prep_filename = ''
        self.frcmod_filename = ''
        self.pdb_filename = ''
        self.original_pdb = ''
        self.map_names_to_old_pdb = {}

        return


def prepare_mols(initial_dir, final_dir, initial_pdb_filename,
                 final_pdb_filename, output_dir, amino=False):
    """

    Args:
        initial_dir (str): Directory containing topology information for
                           initial molecule
        final_dir (str): Directory containing topology information for
                         final molecule
        initial_pdb_filename (str): PDB file for initial molecule
        final_pdb_filename (str): PDB file for final molecule
        output_dir (str): Directory to save output files
        amino (bool): amino?

    Returns:
        MolInfo: Topology, structure and filename info for the initial molecule
        MolInfo: Topology, structure and filename info for the final molecule
    """

    # Get file names of 3 parameter files needed to describe initial molecule
    (initial_prep_filename,
     initial_frcmod_filename,
     initial_ac_filename) = get_param_files(initial_dir)

    # Edit naming for initial molecule to be correctly readable by RDKit,
    # get appropriate name mappings and update ac and prep files accordingly
    # Info on topology (RDKit) and structure (SasMol) stored in MolInfo object
    (initial_mol_info,
     element_counter) = prepare_mol_for_matching(initial_prep_filename,
                                                 initial_ac_filename,
                                                 initial_pdb_filename,
                                                 output_dir,
                                                 amino=amino)

    # Copy frcmod file to be saved with same basename as new ac and prep files
    frcmod_name = initial_mol_info.output_basename + '.frcmod'
    initial_mol_info.frcmod_filename = os.path.join(output_dir, frcmod_name)
    shutil.copyfile(initial_frcmod_filename, initial_mol_info.frcmod_filename)

    # Get the filenames of 3 parameter files needed to describe final molecule
    (final_prep_filename,
     final_frcmod_filename,
     final_ac_filename) = get_param_files(final_dir)

    # Edit naming as for initial molecule
    # Counter ensures atom names unique across both
    (final_mol_info,
     element_counter) = prepare_mol_for_matching(final_prep_filename,
                                                 final_ac_filename,
                                                 final_pdb_filename,
                                                 output_dir,
                                                 counter=element_counter,
                                                 prefix='final_',
                                                 amino=amino)

    # Copy frcmod file to be saved with same basename as edited ac and prep
    # files
    frcmod_name = final_mol_info.output_basename + '.frcmod'
    final_mol_info.frcmod_filename = os.path.join(output_dir, frcmod_name)
    shutil.copyfile(final_frcmod_filename, final_mol_info.frcmod_filename)

    return initial_mol_info, final_mol_info


def get_param_files(param_path):
    """
    Find prep, frcmod and ac files in path provided. If multiple files found
    then give user the choice of which to use. These files provide the Amber
    force field description for a given molecule.

    Args:
        param_path (str): Path to search for the PDB containing atomic
                          coordinates and names

    Returns:
        str: Path to prep file
        str: Path to frcmod file
        str: Path to ac file

    """

    preps = glob.glob(os.path.join(param_path, '*.prep[ic]'))
    preps += glob.glob(os.path.join(param_path, '*.prep'))

    if not preps:

        raise Exception(
            "No ac file found in provided path: {0:s}".format(param_path))

    elif len(preps) == 1:

        prep_filename = preps[0]

    else:

        option_no = get_user_selection("Amber library (prep) file", preps)
        prep_filename = preps[option_no]

    frcmods = glob.glob(os.path.join(param_path, '*.frcmod'))

    if not frcmods:

        raise Exception(
            "Unable to find an ac file in provided path: {0:s}".format(param_path))

    elif len(frcmods) == 1:

        frcmod_filename = frcmods[0]

    else:

        option_no = get_user_selection("Amber frcmod file", frcmods)
        frcmod_filename = frcmods[option_no]

    acs = glob.glob(os.path.join(param_path, '*.ac'))

    if not acs:
        acs = glob.glob(os.path.join(param_path, '*.AC'))

    if not acs:

        raise Exception(
            "No ac file found in provided path: {0:s}".format(param_path))

    elif 'resp.ac' in acs:

        # Default choice for ac file to use
        ac_filename = 'resp.ac'

    elif len(acs) == 1:

        ac_filename = acs[0]

    else:

        option_no = get_user_selection("ac (resp charge) file", acs)
        ac_filename = acs[option_no]

    return prep_filename, frcmod_filename, ac_filename


def prepare_mol_for_matching(prep_filename, ac_filename, pdb_filename,
                             output_dir, counter=None, prefix='init_',
                             amino=False):
    """
    RDKit makes particular requirements of PDBs to be read. Edit input files
    to accommodate this, output edited PDB and read this in to provide
    topology and structure information.

    Args:
        prep_filename (str): Input Amber library filename
        ac_filename (str): Input Amber ac filename (contains resp charge
                           information)
        pdb_filename (str): Input PDB filename
        output_dir (str): Path to output directory
        counter (collections.Counter): Number of each element encountered in
                                       previous molecule(s)
        prefix (str): Text to prepend to output filenames (.pdb, .prep and
                      resp .ac)

    Returns:
        MolInfo: Topology, structure and filename information for the input
                 molecule
        collections.Counter: Updated number of each element encountered

    """

    pdb_basename = os.path.basename(pdb_filename)

    output_basename = prefix + os.path.splitext(pdb_basename)[0]

    output_pdb = os.path.join(output_dir, output_basename + '.pdb')
    output_ac = os.path.join(output_dir, output_basename + '.ac')
    output_prep = os.path.join(output_dir, output_basename + '.prep')

    # Tidy PDB and rename atoms to get a consistent structure which can be
    # read by RDkit (topology) and SasMol (structural information)
    # Obtain objects of both type, mapping of new atom names to originals
    # and count of the element in the ligand
    (rdkit_mol,
     safe_struct,
     name_map,
     element_counter,
     atom_types) = prepare_structure_for_matching(prep_filename,
                                                  pdb_filename,
                                                  output_pdb,
                                                  counter=counter,
                                                  amino=amino)

    cwd = os.getcwd()

    os.chdir(os.path.join(output_dir, 'tmp'))

    # Rename atoms in resp and library files to agree with updated structure
    prepare_param_for_matching(ac_filename, output_pdb, output_ac, output_prep)

    os.chdir(cwd)

    # Combine read information in MolInfo for convenience
    mol_info = MolInfo()

    mol_info.original_pdb = pdb_filename

    mol_info.pdb_filename = output_pdb
    mol_info.ac_filename = output_ac
    mol_info.prep_filename = output_prep

    mol_info.mol = rdkit_mol
    mol_info.struct = safe_struct

    mol_info.atom_types = atom_types

    mol_info.map_names_to_old_pdb = name_map
    mol_info.output_basename = output_basename

    return mol_info, element_counter


def prepare_param_for_matching(ac_filename, naming_pdb_filename,
                               output_ac, output_prep):
    """
    Update the input ac and prep files (containing resp charge and topology
    information respectively) to us the naming conventions in
    naming_pdb_filename. Output updated files.

    Args:
        ac_filename (str): Path to the input ac file
        naming_pdb_filename (str): Path to PDB file with updated atom naming
        output_ac (str): Path to use for output ac file
        output_prep (str): Path to use for output Amber library file

    Returns:

    """

    returncode = subprocess.call(['antechamber', '-fi', 'ac',
                                  '-fo', 'ac', '-i', ac_filename,
                                  '-o', output_ac,
                                  '-ao', 'name',
                                  '-a', naming_pdb_filename,
                                  '-fa', 'pdb'])

    if returncode:
        print(
            'Antechamber error: Unable to convert {0:s} to naming from {1:s}'.format(
                ac_filename,
                naming_pdb_filename))
        sys.exit(1)

    returncode = subprocess.call(['antechamber', '-i', output_ac,
                                  '-fi', 'ac',
                                  '-o', output_prep,
                                  '-fo', 'prepi'], shell=True)

    if returncode:
        print(
            'Antechamber error: Unable to convert {0:s} to prep file'.format(output_ac))
        sys.exit(1)

    return


def prepare_structure_for_matching(prep_filename, pdb_filename,
                                   output_pdb, counter=None, amino=False):
    """
    Create an RdKit Mol object from the input structure. The atoms of the
    input structure are renamed to ensure they can be read by RdKit. The
    RdKit Mol, edited structure and mapping between atom names are
    returned. A PDB file using the new naming is saved to file.

    Args:
        prep_filename (str): Path to the library (prep) file
        pdb_filename (str): Path to the PDB containing atomic coordinates
                            and names
        output_pdb (str): Path to the output (atom renamed) PDB file

    Kwargs:
        counter (collections.Counter): Count of each element type seen. Allows
                                       atoms to be uniquely identified by name
                                       across molecules

    Returns:
        rdkit.Chem.rdchem.Mol: Contains molecule atoms and connections
        sasmol.SasMol: Structure containing renamed atoms
        dict: Map between input and output atom names
        collections.Counter: Updated element counts
        dict: Index to atom type mapping

    """

    prep, original_structure = read_pdb_prep_pair(prep_filename, pdb_filename)

    atom_names = original_structure.name()

    atom_types = {}
    for idx in range(original_structure.natoms()):
        atom_name = atom_names[idx]
        try:
            atom_types[idx] = prep.atom_type_from_name[atom_name]
        except:
            test_name = atom_name[1:] + atom_name[0]
            atom_types[idx] = prep.atom_type_from_name[test_name]

    if counter:

        (safe_struct,
         name_map,
         element_counter) = create_safe_atom_name_pdb(original_structure,
                                                      output_pdb,
                                                      counter=counter,
                                                      amino=amino)

    else:

        (safe_struct,
         name_map,
         element_counter) = create_safe_atom_name_pdb(original_structure,
                                                      output_pdb,
                                                      amino=amino)

    rdkit_mol = Chem.MolFromPDBFile(output_pdb, removeHs=False)

    return rdkit_mol, safe_struct, name_map, element_counter, atom_types


def create_safe_atom_name_pdb(structure, output_pdb_name, amino=False,
                              **kwargs):
    """
    Rename the atoms in the input structure so that they are safe for Amber
    and RdKit usage. Output to a new PDB.

    Args:
        structure (sasmol.SasMol): Contains information from original PDB
        output_pdb_name (str): Filename for output PDB with edited atom names

    Returns:
        sasmol.SasMol: Contains information for renamed PDB
        dictionary: Map new atom names to the original ones

    """

    # Allow element counter initialization so atoms in multiple structures all
    # have unique ids
    if 'counter' in kwargs:
        element_count = kwargs['counter']
    else:
        element_count = Counter()

    new_names = []
    new_elements = []

    original_names = structure.name()

    # Copy molecule to new structure
    new_structure = copy.deepcopy(structure)

    # Rename atoms as XN where X = element and N = no. element counted
    for idx in range(structure.natoms()):
        # HACK: Derive element from atom name
        # SasMol doing this already but confused by multi-letter elements
        # This is a temporary correction
        element = element_from_pdb_atom_name(original_names[idx], amino=amino)

        element_count.update([element])

        new_name = '{0:s}{1:d}'.format(element, element_count[element])

        new_names.append(new_name)
        new_elements.append(element)

    # Update structure naming and save as PDB
    new_structure.setName(new_names)
    new_structure.setElement(new_elements)

    new_structure.write_pdb(output_pdb_name, 0, 'w')

    name_map = OrderedDict(zip(new_names, original_names))

    return new_structure, name_map, element_count


def element_from_pdb_atom_name(atom_name, amino=False):
    """
    Try to extract the element name from the atom name provided (from a PDB)

    Args:
        atom_name (str): Atom name read from a PDB file

    Returns:
        str: Element name (in upper case)

    """

    name_element = atom_name.strip('0123456789').upper()

    # Standard elements
    if name_element in ['C', 'O', 'N', 'H', 'P', 'K', 'S', 'SI', 'BR', 'CL',
                        'FL', 'F', 'NA', 'ZN', 'FE']:

        element = name_element

    # So if we were sealing with protein like chains CA would be a Carbon alpha
    # Mostly will not be so, assume CA is what is actually meant
    elif name_element == 'CA':

        if amino:
            element = 'C'
        else:
            print('WARNING: CA atom name interpretted as Calcium\n')
            element = name_element

    elif name_element == 'HN':

        print('WARNING: HN atom name interpretted as Hydrogen\n')
        element = 'H'

    # Terminal oxygens
    elif name_element in ['OXT', 'OT']:

        element = 'O'

    # For protein like Carbons
    elif name_element in ['CB', 'CD', 'CG']:

        element = 'C'

    else:
        raise Exception(
            'Unable to determine element for atom {0:s}'.format(atom_name))

    return element

