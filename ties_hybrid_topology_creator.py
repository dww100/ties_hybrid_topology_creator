#! /usr/bin/env python
"""
Build a hybrid Amber topology from single topologies for two ligands ready for
use in TI calculations.
Elements of the initial molecule 'disappear' whilst those of the final 'appear'
"""
from __future__ import print_function

import sys, os
import subprocess
import argparse
import shutil
import glob
import itertools
import copy
import numpy as np

from collections import OrderedDict
from collections import Counter

from rdkit import Chem
from rdkit.Chem import rdFMCS

import sasmol.sasmol as sasmol

from bac_ties.bac_amber_utils import *
from bac_ties.bac_atom_data import *

input = raw_input

class MolInfo():
    """
    Collects molecular information and associated filenames.

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
        raise Exception('Unable to determine element for atom {0:s}'.format(atom_name))

    return element


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

    Returns:
        MolInfo: Topology, structure and filename info for the initial molecule
        MolInfo: Topology, structure and filename info for the final molecule
    """

    # Get filenames of 3 parameter files needed to describe initial molecule
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

    # Copy frcmod file to be saved with same basename as edited ac and prep files
    frcmod_name = final_mol_info.output_basename + '.frcmod'
    final_mol_info.frcmod_filename = os.path.join(output_dir, frcmod_name)
    shutil.copyfile(final_frcmod_filename, final_mol_info.frcmod_filename)

    return initial_mol_info, final_mol_info


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

    output_pdb = os.path.join(output_dir, output_basename+ '.pdb')
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
        print('Antechamber error: Unable to convert {0:s} to naming from {1:s}'.format(ac_filename, naming_pdb_filename))
        sys.exit(1)

    returncode = subprocess.call(['antechamber', '-i', output_ac,
                              '-fi', 'ac',
                              '-o', output_prep,
                              '-fo', 'prepi'])

    if returncode:
        print('Antechamber error: Unable to convert {0:s} to prep file'.format(output_ac))
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
        atom_types[idx] = prep.atom_type_from_name[atom_name]

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


def get_user_selection(text, option_list, header='', default=0):
    """
    Allow user to select an option from list of options.

    Note: this is only working for string options at present

    Args:
        text (str): Text describing what is being chosen
        option_list (list): List of options from which to select

    Kwargs:
        header (str): Description for columned input
        default (int): Default option text (-1 = no default)

    Returns:
        int: index of selected option in list
    """

    print("Select {0:s}:".format(text))

    print(header)
    for item in enumerate(option_list):
        print('[{0:d}] {1:s}'.format(item[0], item[1]))

    selected = False

    while not selected:

        try:

            if default != -1:
                idx = int(input("Choice [Standard = {0:d}]: ".format(default)))
            else:
                idx = int(input("Choice: "))

            if 0 <= idx < len(option_list):

                selected = True

            else:

                print("You must select one of the presented options")

        except ValueError:

            print("Invalid input")

    return idx


def get_user_submatch(submatches, atom_info, q_diffs, q_max_atom_diffs):
    """
    Get user to select the region to be used as common atoms in the hybrid
    ligand from the selected submatches.

    Args:
        submatches (list): List of lists of atom indices in each submatch
        atom_info (dict): Naming and connectivity (AtomInfo) information
                          by atom index
        q_diffs (list): Difference in charge between submatch in initial and
                        final ligands
        q_max_atom_diffs (list): Maximum individual atom charge difference
                                (initial to final ligand) for each submatch

    Returns:
        int: Index of selected submatch in input list

    """

    option_list = []

    for idx in range(len(submatches)):

        option = ','.join(atom_info[x].name for x in submatches[idx])

        option_list.append('{0:f}\t{1:f}\t{2:s}'.format(q_diffs[idx],
                                                        q_max_atom_diffs[idx],
                                                        option))

    header_txt = "#\tq diff\tatom q diff\tatom list"

    selected = get_user_selection('substructure match',
                                  option_list,
                                  header=header_txt,
                                  default=-1)

    return selected


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


def remove_ring_error(mol, matched, atom_info=None):
    """
    Remove incomplete rings from lists of matched atoms.

    Args:
        mol (rdkit.Chem.rdchem.Mol): Contains information on atoms in molecule
        matched (list): Matched atom indices
        atom_info (dict): Naming and connectivity (AtomInfo) information
                          by atom index

    Returns:
        list: list of matches atom indices with partial ring selections
              removed

    """

    atom_rings = mol.GetRingInfo().AtomRings()


    all_ring_atoms = [x for sublist in atom_rings for x in sublist]

    updated_match = []

    # Include only atoms matching complete rings in updated list
    for ring in atom_rings:
        if set(ring).issubset(matched):
            updated_match += ring

    # Add atoms not in rings at all
    updated_match += list(set(matched).difference(all_ring_atoms))

    if atom_info:

        # In cases where ring links multiple matching sections
        # need to find largest connected section

        disjoint_sections = []

        for idx in updated_match:

            if not disjoint_sections:

                disjoint_sections.append([idx])

            else:

                # Add atom to either a new list or a list containing an atom
                # it is bound to
                info = atom_info[idx]

                section_idx = None

                bound_match = list(set(info.bound).intersection(set(updated_match)))

                for bound_idx in bound_match:

                    for check_idx in range(len(disjoint_sections)):

                        if bound_idx in disjoint_sections[check_idx]:

                            section_idx = check_idx
                            break

                # Add this atom and any bound to it to appropriate list
                if section_idx == None:

                    tmp = [idx]
                    tmp += bound_match
                    disjoint_sections.append(tmp)

                else:

                    disjoint_sections[section_idx].append(idx)
                    disjoint_sections[section_idx] += bound_match

        # Due to atom/bond ordering some linked sections may be
        # separated - combine these and retain only unique indices
        final_disjoint_sections = []
        while len(disjoint_sections) > 0:

            first, rest = disjoint_sections[0], disjoint_sections[1:]
            first = set(first)

            len_previous_first = -1

            while len(first) > len_previous_first:

                len_previous_first = len(first)

                rest2 = []

                for existing_list in rest:

                    if len(first.intersection(set(existing_list))) > 0:

                        first |= set(existing_list)

                    else:

                        rest2.append(existing_list)

                rest = rest2

            final_disjoint_sections.append(first)
            disjoint_sections = rest

        # Select the largest disjoint section
        updated_match = max(final_disjoint_sections, key=len)

    return updated_match


def get_bridge_atoms(mol, matched_idxs):
    """
    Produce list of atoms that bridge between matched region and the rest of
    the molecule.

    Args:
        mol (rdkit.Chem.rdchem.Mol): Molecular connectivity information
        matched_idxs (list): Indices of matched atoms

    Returns:
        list: Indices for matched atoms bound to non-matched atoms

    """

    bridge_atoms = []

    for bond in mol.GetBonds():

        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()

        begin_matched = begin_idx in matched_idxs
        end_matched = end_idx in matched_idxs

        if begin_matched != end_matched:

            if begin_matched and begin_idx not in bridge_atoms:

                bridge_atoms.append(begin_idx)

            elif end_idx not in bridge_atoms:

                bridge_atoms.append(end_idx)


    return bridge_atoms


def get_stop_ring_remove(idx, atom_info, atom_rings, matched):
    """
    Get indices of non-ring atoms connected to the ring containing the atom
    with index idx (or rings connected to that ring). These atoms are where
    removal of a ring from the matched region stops. Due to this goal only
    indices of matched atoms are reported.

    Args:
        idx (int): Atom index for atom being removed
        atom_info (dict): Naming and connectivity (AtomInfo) information
                          by atom indexm
        atom_rings (list): List of lists of atom indices for all rings in
                           molecule
        matched (list): Atom indices for the macthed region of the molecule

    Returns:
        list: indices of atoms where ring removal finishes

    """

    stop_idxs = []

    # Start with the original atom being removed
    start_info = atom_info[idx]

    # Get indices of all rings of which it is part
    for ring_no in start_info.rings:

        for atom_idx in atom_rings[ring_no]:

            info = atom_info[atom_idx]

            for bound_idx in info.bound:

                if not atom_info[bound_idx].rings:

                    if (bound_idx not in stop_idxs) and (bound_idx in matched):

                        stop_idxs.append(bound_idx)

    return stop_idxs


def get_linked_to_remove(idx, atom_info, matched, stop_list):
    """
    Obtain list of atoms to be removed from the macthed list in order to
    cleanly remove atom with index idx. This includes atoms in rings
    containing idx (and rings linked directly to these) and atoms bound
    to atom idx and nothing else. The process is recursive to allow the
    removal of sections of branched structures.

    Args:
        idx (int): Index of atom to be removed from match
        atom_info (dict): Naming and connectivity (AtomInfo) information
                          by atom index
        matched (list): Indices of matched atoms
        stop_list (list): Indices of atoms that should end removal region.

    Returns:
        list: Indices to remove from match

    """

    to_remove = [idx]

    stop_list.append(idx)

    atom = atom_info[idx]

    for bonded_idx in atom.bound:

        if (bonded_idx in matched) and (bonded_idx not in stop_list):

            to_remove += get_linked_to_remove(bonded_idx, atom_info,
                                              matched, stop_list)

    return to_remove


def get_submatches(mol, atom_info, match_idxs):
    """
    Obtain a list of submatches within the provided selection of
    matching atoms within mol. The removal of any atom in a ring
    from the match results in the removal of the entire ring.

    Args:
        mol (rdkit.Chem.rdchem.Mol): Molecular connectivity information
        atom_info (dict): Naming and connectivity (AtomInfo) information
                          by atom index
        match_idxs (list): Inidices of matched atoms

    Returns:
        list: List of lists of atom indices for each submatch identified

    """

    atom_rings = mol.GetRingInfo().AtomRings()

    options = []
    trial_bridges = []

    # Identify atoms which link matched to unmatched region of the molecule
    original_bridge_idxs = get_bridge_atoms(mol, match_idxs)

    # Starting from each bridge between matched and unmatched region list
    # segments generated by removing each atom sequentially from the match
    # Removal of an atom in a ring removes the ring and attached atoms
    # are removed alongside the atom of interest
    for start_bridge_idx in range(len(original_bridge_idxs)):

        # List to store potential sections to starting from each bridge in
        # the original match
        options.append([])

        # Start search at bridge atom from the full match
        bridge_idxs = [original_bridge_idxs[start_bridge_idx]]

        # Loop until don't find any new bridge indices
        while True:

            for bridge_idx in bridge_idxs:

                bridge_info = atom_info[bridge_idx]

                # Get list of atoms where region removal should stop
                if bridge_info.rings:

                    potential_stop_idxs = get_stop_ring_remove(bridge_idx,
                                                               atom_info,
                                                               atom_rings,
                                                               match_idxs)

                else:

                    potential_stop_idxs = [x for x in bridge_info.bound if x in match_idxs]

                # Get selections of atoms to remove for each potential section
                # to be removed
                for stop_idx in potential_stop_idxs:

                    stop_idxs = [stop_idx]

                    remove_selection = get_linked_to_remove(bridge_idx,
                                                            atom_info,
                                                            match_idxs,
                                                            stop_idxs)

                    remove_selection.sort()

                    # Add unique options to the list for starting bridge atom
                    if ((len(remove_selection) < len(match_idxs) - 1) and
                        (remove_selection not in options[start_bridge_idx])):

                        options[start_bridge_idx].append(remove_selection)

            # Record which bridge indices have already been evaluated
            trial_bridges += bridge_idxs

            # Reset list of bridges to be evaluated and reset to contain those
            # created when each section identified above is removed
            bridge_idxs = []

            for option in options[start_bridge_idx]:

                new_match = list(set(match_idxs) - set(option))

                # Don't re-examine bridges that we have checked already
                bridge_idxs += [x for x in get_bridge_atoms(mol, new_match) if x not in trial_bridges]

            if len(set(bridge_idxs)) < 1:
                break

    # Need a comprehensive list of all possible removals that can be used to
    # produce a valid submatch
    # a) Create all possible combinations of section removals
    tmp = []

    for comb_length in range(1, len(options)+1):

        for combination in itertools.combinations(options, comb_length):

            for product in itertools.product(*combination):

                combined_removal_atoms = set([item for sublist in product for item in sublist])

                if combined_removal_atoms not in tmp:

                    tmp.append(combined_removal_atoms)

    # Get a unique set of options (sorted to make comparisons easier later on)
    full_options = [sorted(list(x)) for x in set(tuple(x) for x in tmp)]

    # b) Ensure single bridge atoms are recorded
    for bridge_options in options:

        for option in bridge_options:

            option.sort()

            if option not in full_options:

                full_options.append(option)

    # c) Add full match to options (i.e. removal of no atoms from match)
    full_options.append([])

    full_options.sort(key=len)

    # Submatches = match - removed atoms for each option created above
    submatches = []

    for removal_list in full_options:

        submatch = list(set(match_idxs) - set(removal_list))

        submatch.sort()

        if submatch:
            submatches.append(submatch)

    # for i in range(len(submatches)):
    #     print(submatches[i])
    #     print(full_options[i])

    # sys.exit()

    return submatches


def get_charge_diffs_submatches(submatches, matched_idx_map,
                                initial_charge_map, final_charge_map):
    """
    Get overal land maximum individual atom charge differences between the
    initial and final molecules for each submatch passed in.

    Args:
        submatches (list): List of lists of atom indices for each submatch
        matched_idx_map (dict): Map index of atom in initial to that in the
                                final molecule
        initial_charge_map (bac_ties.bac_atom_data.PdbRdkitChargeMap): Charges for atom indices for initial molecule
        final_charge_map (bac_ties.bac_atom_data.PdbRdkitChargeMap): Charges for atom indices for final molecule

    Returns:
        list: Overall charge differences for each submatch
        list: Maximum atomic charge difference for each submatch

    """

    q_diffs = []
    q_max_atom_diffs = []

    for submatch in submatches:

        q_diff_total = 0
        q_max_atom_diff = 0

        for idx in submatch:

            q_init = initial_charge_map.charge_rdkit_idx(idx)

            idx_final = matched_idx_map[idx]

            q_final = final_charge_map.charge_rdkit_idx(idx_final)

            q_diff = q_init - q_final

            if abs(q_diff) > q_max_atom_diff:
                q_max_atom_diff = abs(q_diff)

            q_diff_total += q_diff

        q_max_atom_diffs.append(q_max_atom_diff)
        q_diffs.append(abs(q_diff_total))

    return q_diffs, q_max_atom_diffs


def output_submatches_file(submatches, selected,
                           q_diffs, q_max_atom_diffs,
                           initial_info, final_info,
                           matched_idx_map, output_dir):
    """
    Output a summary file detailing the possible submatches and the atomic
    charge differences between the initial and final molecules

    Args:
        submatches (list): List of lists of atom indices for each submatch
        selected (int): Index of selected submatch in teh list above
        q_diffs (list): Overall charge difference for each submatch
        q_max_atom_diffs (list): Maximum individual atom charge difference
                                 for each submatch
        initial_info (dict): Naming and connectivity information (AtomInfo)
                             by atom index (initial molecule)
        final_info (dict): Naming and connectivity information (AtomInfo)
                           by atom index (final molecule)
        matched_idx_map (dict): Map index of atom in initial to that in the
                                final molecule
        output_dir (str): Path in which to save output

    Returns:

    """

    out_filename = os.path.join(output_dir, 'original_matches.txt')

    out_file = open(out_filename, 'w')

    print('List of all available matches ', file=out_file)
    print('Uses atom names from initial molecule', file=out_file)
    print('Atom_Q_diff = maximum difference in selection\n', file=out_file)

    print('#\tQ_diff\tAtom_Q_diff\tMatched atoms', file=out_file)

    for sub_no in range(len(submatches)):

        match = submatches[sub_no]
        option = ' '.join(initial_info[x].name for x in match)

        print('{0:d}\t{1:f}\t{2:f}\t{3:s}'.format(sub_no, q_diffs[sub_no],
                                                  q_max_atom_diffs[sub_no],
                                                  option),
              file=out_file)

    print(file=out_file)

    print('Selected: {0:d}'.format(selected), file=out_file)

    print(file=out_file)

    print('Name_initial\tName_final\tQ_initial\tQ_final\tAtom_Q_diff', file=out_file)

    for idx1, idx2 in matched_idx_map.iteritems():

        name1 = initial_info[idx1].name
        name2 = final_info[idx2].name
        charge1 = initial_info[idx1].charge
        charge2 = final_info[idx2].charge
        diff = charge1 - charge2

        print('{0:s}\t{1:s}\t{2:f}\t{3:f}\t{4:f}'.format(name1, name2,
                                                         charge1, charge2,
                                                         diff),
              file=out_file)

    print('\nNote: Names from initial are used in final output', file=out_file)

    print('\nUnmatched atoms from initial:', file=out_file)

    unmatched_idxs = [initial_info[x].name for x in initial_info.keys() if x not in matched_idx_map.keys()]
    unmatched_txt = ' '.join(unmatched_idxs)
    print(unmatched_txt, file=out_file)

    print('\nUnmatched atoms from final:', file=out_file)

    unmatched_idxs = [final_info[x].name for x in final_info.keys() if x not in matched_idx_map.keys()]
    unmatched_txt = ' '.join(unmatched_idxs)
    print(unmatched_txt, file=out_file)

    out_file.close()

    return


def check_atom_type_match(matched_idx_map, initial_mol_info, final_mol_info):
    """
    Check atom types for matched atoms are the same in both molecules. If not
    the same report to user and ask if they wish to continue and ends if not.

    Args:
        matched_idx_map (dict): Map index for initial molecule to index of the
                                corresponding atom in the final molecule
        initial_mol_info (MolInfo()): Structure, name and type information for
                                      the initial molecule
        final_mol_info (MolInfo()): Structure, name and type information for
                                    the final molecule

    Returns:

    """

    init_names = initial_mol_info.struct.name()
    final_names = final_mol_info.struct.name()

    type_match_error = False

    for init_idx, final_idx in matched_idx_map.iteritems():

        init_type = initial_mol_info.atom_types[init_idx]
        final_type = final_mol_info.atom_types[final_idx]

        if init_type != final_type:
            init_name = init_names[init_idx]
            final_name = final_names[final_idx]

            print('WARNING: Matched atoms of different types: {0:s} ({1:s}) - {2:s} ({3:s})'.format(init_name,
                                                                                                    init_type,
                                                                                                    final_name,
                                                                                                    final_type))
            type_match_error = True

    if type_match_error:

        choice = input("Given type match error do you want to proceed (y/[n])?")
        if choice.lower() != 'y':
            sys.exit(1)

    return


def compare_ligands(initial_dir, initial_pdb, final_dir, final_pdb,
                    output_dir, manual, tolerance=0.1, atom_tolerance=0.1,
                    amino=False):
    """
    Compare ligands and select a matching region from both that will be
    used as the common region when building the final hybrid ligand.
    Return information on both ligands and the matched region.

    Args:
        initial_dir (str): Path to files describing initial molecule
        initial_pdb (str): Path to PDB for initial molecule
        final_dir (str): Path to files describing final molecule
        final_pdb (str): Path to PDB for final molecule
        output_dir (str): Path to save output
        manual (bool): Will teh user manually select the submatch to use?

    Returns:
        MolInfo: Topology, structure and filename information for the
                 initial molecule
        dict: Atomic information for initial molecule: keys = atom index
              (RDKit style), values = AtomInfo containing bond and naming
              information
        MolInfo: Topology, structure and filename information for the
                 final molecule
        dict: Atomic information for final molecule: keys = atom index
              (RDKit style), values = AtomInfo containing bond and naming
              information
        dict: Map of RDKit indices from initial to final molecule for
              matched atoms
        list: Indices of atoms of the selected submatch (uses initial molecule
              and RDKit style indexing.

    """

    # Read structure and topology information from PDBs
    # Identify matching ac, prep and frcmod files in input directories
    # RDKit requires atom name changes, update input files & move to output_dir
    initial_mol_info, final_mol_info = prepare_mols(initial_dir, final_dir,
                                                    initial_pdb, final_pdb,
                                                    output_dir, amino=amino)

    # Map charges from resp ac file to atom names and indexes
    initial_charge_map = create_charge_idx_map(initial_mol_info.ac_filename,
                                               initial_mol_info.struct)

    final_charge_map = create_charge_idx_map(final_mol_info.ac_filename,
                                             final_mol_info.struct)

    # Get dictionary containing AtomInfo for each atom
    # (uses RDKit index as key)
    initial_atom_info = get_atom_info_rdkit(initial_mol_info.mol,
                                            charge_map=initial_charge_map)

    final_atom_info = get_atom_info_rdkit(final_mol_info.mol,
                                          charge_map=final_charge_map)

    # Find largest matching region between PDBs (Miximal Common Substructure)

    rdkit_mols = [initial_mol_info.mol, final_mol_info.mol]
    mcs_result = rdFMCS.FindMCS(rdkit_mols, completeRingsOnly=True)

    if mcs_result:

        shared = Chem.MolFromSmarts(mcs_result.smartsString)

        initial_match_idxs = initial_mol_info.mol.GetSubstructMatch(shared)
        final_match_idxs = final_mol_info.mol.GetSubstructMatch(shared)

        # Map to get equivalent atom index in final molecule to the in the
        # initial model (matched region)
        matched_idx_map = dict(zip(initial_match_idxs, final_match_idxs))

        initial_match_idxs = remove_ring_error(initial_mol_info.mol,
                                               initial_match_idxs,
                                               atom_info=initial_atom_info)

        # Update final matched idx list and initial to final map to account
        # for removed ring atoms

        final_match_idxs = [matched_idx_map[x] for x in initial_match_idxs]

        matched_idx_map = dict(zip(initial_match_idxs, final_match_idxs))

        if len(initial_match_idxs) == 0:
            print("No matching region (excluding incomplete rings) could be identified!")
            sys.exit(0)




        # Get all potential sub-matches with complete rings
        # (includes complete match)
        submatches = get_submatches(initial_mol_info.mol, initial_atom_info,
                                    initial_match_idxs)

        (q_diffs,
         q_max_atom_diffs) = get_charge_diffs_submatches(submatches,
                                                         matched_idx_map,
                                                         initial_charge_map,
                                                         final_charge_map)

        # Manual selection of submatch by user or automatic choice
        if manual:

            selected_idx = get_user_submatch(submatches, initial_atom_info,
                                             q_diffs, q_max_atom_diffs)

        else:


            all_options = range(len(submatches))
            filtered_q_total = [x for x in all_options if q_diffs[x] < tolerance]
            filtered_q_atom = [x for x in filtered_q_total if q_max_atom_diffs[x] < atom_tolerance]

            all_filtered_matches = []

            if len(filtered_q_atom) > 0:

                for submatch_idx in filtered_q_atom:

                    submatch = submatches[submatch_idx]

                    no_heavy = 0

                    for idx in submatch:

                        if initial_atom_info[idx].name.strip('0123456789') != 'H':
                            no_heavy += 1

                    if no_heavy > 1:
                        all_filtered_matches.append(submatch_idx)

                if len(all_filtered_matches) > 0:

                    selected_idx = all_filtered_matches[0]

                else:

                    print('Warning: Only acceptable matches with at most 1 heavy atom found')
                    selected_idx = filtered_q_atom[0]

            else:

                print("No options fit standard filter - select manually\n")

                selected_idx = get_user_submatch(submatches, initial_atom_info,
                                                 q_diffs, q_max_atom_diffs)

    else:

        print("No matching region could be identified!")
        sys.exit(1)

    selected_submatch = submatches[selected_idx]

    # Output all submatches and information informing choice
    output_submatches_file(submatches, selected_idx,
                           q_diffs, q_max_atom_diffs,
                           initial_atom_info, final_atom_info,
                           matched_idx_map, output_dir)

    return (initial_mol_info, initial_atom_info,
            final_mol_info, final_atom_info, matched_idx_map, selected_submatch)


def rename_common_atoms_final(final_struct, initial_submatch_idxs,
                              initial_atom_info, idx_map):
    """
    Rename atoms in the input final structure to reflect those for the
    equivalent atoms in the initial molecule.

    Args:
        final_struct (sasmol.SasMol): Structure information for final molecule
        initial_submatch_idxs (list): Indices of chosen matching atoms in
                                      initial molecule
        initial_atom_info (dict): Naming and connectivity (AtomInfo)
                                  information by atom index
        idx_map (dict): Mapping of initial to final indices for matched atoms

    Returns:

    """

    final_common_idxs = [idx_map[x] for x in initial_submatch_idxs]
    inverse_idx_map = {v: k for k, v in idx_map.items()}

    current_names = final_struct.name()
    new_names = []

    for idx in range(final_struct.natoms()):

        if idx in final_common_idxs:

            initial_idx = inverse_idx_map[idx]
            new_name = initial_atom_info[initial_idx].name
            new_names.append(new_name)

        else:

            new_names.append(current_names[idx])


    final_struct.setName(new_names)

    return


def update_final_description(final_mol_info, initial_submatch_idxs,
                             initial_atom_info, idx_map, output_dir):
    """
    Update the files (.pdb, .prep and .ac) describing the final ligand to use
    naming consistent with matched area in the initial molecule.

    Args:
        final_mol_info (MolInfo): Topology, structure and filename information
                                  for the final molecule
        initial_submatch_idxs (list): Indices of selected matched atoms in
                                      initial molecule
        initial_atom_info (dict): Atomic charge and naming information by index
        idx_map (dict): Map of initial to final molecule indices
        output_dir (str): Path to store output files

    Returns:

    """

    renamed_pdb = os.path.join(output_dir, 'final_renamed.pdb')
    renamed_ac = os.path.join(output_dir, 'final_renamed.ac')
    renamed_prep = os.path.join(output_dir, 'final_renamed.prep')

    rename_common_atoms_final(final_mol_info.struct,
                              initial_submatch_idxs,
                              initial_atom_info,
                              idx_map)

    final_mol_info.struct.write_pdb(renamed_pdb, 0, 'w')
    final_mol_info.mol = Chem.MolFromPDBFile(renamed_pdb)

    cwd = os.getcwd()

    os.chdir(os.path.join(output_dir, 'tmp'))

    prepare_param_for_matching(final_mol_info.ac_filename,
                               renamed_pdb, renamed_ac, renamed_prep)

    os.chdir(cwd)

    final_mol_info.ac_filename = renamed_ac
    final_mol_info.prep_filename = renamed_prep
    final_mol_info.pdb_filename = renamed_pdb

    return


def calc_coor_com(struct, atom_names, frame=0):
    """
    Calculate the coordinates and centre of mass of the selected atoms in
    struct. Coordinates are returned in the order of the names in the
    input list.

    Args:
        struct (sasmol.SasMol: Molecular structure information
        atom_names (list): Atom names to use in selection from struct
        frame (int): Number of frame to use for coordinates in struct

    Returns:
        np.array: Array containing x,y,z coordinates for each selected atoms
        np.array: Coordinates of the centre of mass of teh selected atoms

    """


    formatted_atom_names = seg_list = ['"{0:s}"'.format(x) for x in atom_names]
    filter_txt = 'name[i] in [{0:s}]'.format(','.join(formatted_atom_names))

    err, mask = struct.get_subset_mask(filter_txt)
    tmp_struct = sasmol.SasMol(0)
    struct.copy_molecule_using_mask(tmp_struct, mask, frame)

    com = tmp_struct.calccom(frame)

    # Coordinates need to be in the order of the input atom names
    coor = []
    tmp_coor = tmp_struct.coor()[frame]
    tmp_names = tmp_struct.name()

    for name in atom_names:

        idx = tmp_names.index(name)

        coor.append(tmp_coor[idx])

    coor = np.array(coor)

    return coor, com


def align_molecules_using_common(align_struct, target_struct,
                                 common_atom_names, frame=0):
    """
    Align structure on target using the selected common atom names.

    Args:
        align_struct (sasmol.SasMol): Atomic information for structure to be
                                      aligned
        target_struct (sasmol.SasMol): Atomic information for target structure
                                       for alignment
        common_atom_names (list): Atom names to be used in the alignment
        frame (int): Number of frame to use for coordinates in alignment

    Returns:

    """

    align_coor, align_com = calc_coor_com(align_struct, common_atom_names)
    target_coor, target_com = calc_coor_com(target_struct, common_atom_names)

    align_struct.align(frame, align_coor, align_com, target_coor, target_com)

    return


def calculate_average_charges(initial_atom_info, final_atom_info,
                              submatch, matched_idx_map):
    """
    Calculate the average charge of each atom in the submatch, using the
    charges in the initial and final molecules.

    Args:
        initial_atom_info (dict): Atomic charge and naming information by index
                                  (initial molecule)
        final_atom_info (dict): Atomic charge and naming information by index
                                (final molecule)
        submatch (list): Atomic indices of chosen matching atoms
        matched_idx_map: Map of initial to final molecule indices

    Returns:
        list: Average charge of each atom in submatch

    """

    average_charges = {}

    for initial_idx in submatch:

        atom_name = initial_atom_info[initial_idx].name
        final_idx = matched_idx_map[initial_idx]

        charge1 = initial_atom_info[initial_idx].charge
        charge2 = final_atom_info[final_idx].charge

        average_charges[atom_name] = (charge1 + charge2) / 2.0

    return average_charges


def write_charge_constraint_file(avg_charges, struct, atom_info, filename):
    """
    Write a constraint file for use in Antechamber that constrains the charge
    of each atom for which an average charge is provided to this charge.

    Args:
        avg_charges (dict): Average charge by atom name (matched atoms only)
        struct (sasmol.SasMol: Atomic information for molecule
        atom_info (dict): Atomic charge and naming information by index
        filename (str): Filename into which to save the constraints

    Returns:

    """

    out_file = open(filename, 'w')

    for idx, atom_name in enumerate(struct.name()):

        if atom_name in avg_charges:

            avg_charge = avg_charges[atom_name]

            # idx use here is a bit dubious - but should work
            original_charge = atom_info[idx].charge

            charge_diff = abs(avg_charge-original_charge)

            if charge_diff > 0.5:
                print("Charge difference > 0.5 (0:f) for atom {1:s}".format(charge_diff, atom_name))

            print('CHARGE {0:.6f} {1:d} {2:s}'.format(avg_charge,
                                                      idx + 1,
                                                      atom_name),
                  file=out_file)

    out_file.close()

    return


def create_updated_prep_frcmod(src_dir, ac_filename,
                               constraint_filename, mol_name):
    """
    Use antechamber to produce prep and frcmod files for selected molecule,
    constraining selected atomic charges to average values (between initial
    and final molecules).

    Args:
        src_dir (str): Path in which to find ESP file
        ac_filename (str): Path to original ac file
        constraint_filename (str): Path to charge constraint file
        mol_name (str): Name to be given to output molecule

    Returns:

    """

    esp_filename = os.path.join(src_dir, 'ANTECHAMBER.ESP')

    returncode = 0

    # Get input files for two stage resp fitting (using constraints)
    returncode += subprocess.call(['respgen', '-i', ac_filename,
                                  '-o', mol_name + '.respin1',
                                  '-f', 'resp1',
                                  '-e', '1',
                                  '-a', constraint_filename])


    returncode += subprocess.call(['respgen', '-i', ac_filename,
                                  '-o', mol_name + '.respin2',
                                  '-f', 'resp2',
                                  '-e', '1',
                                  '-a', constraint_filename])

    qin_filename = mol_name + '_QIN'
    os.rename('QIN', qin_filename)

    # Perform resp fitting
    returncode += subprocess.call(['resp', '-O', '-i', mol_name + '.respin1',
                                  '-o', mol_name + '.respout1',
                                  '-e', esp_filename,
                                  '-t', mol_name + '_qout',
                                  '-q', qin_filename])

    returncode += subprocess.call(['resp', '-O', '-i', mol_name + '.respin2',
                                  '-o', mol_name + '.respout2',
                                  '-e', esp_filename,
                                  '-t', mol_name + '_QOUT',
                                  '-q', mol_name + '_qout'])

    # Convert resp output into ac file format
    returncode += subprocess.call(['antechamber', '-i', ac_filename,
                                  '-fi', 'ac', '-c', 'rc',
                                  '-cf', mol_name + '_QOUT',
                                  '-o', mol_name + '.prep',
                                  '-fo', 'prepi',
                                  '-rn', mol_name.upper()])

    # Output missing force field parameters in frcmod format
    returncode += subprocess.call(['parmchk', '-i', mol_name + '.prep',
                                  '-f', 'prepi',
                                  '-o', mol_name + '.frcmod'])

    if returncode:

        print('ERROR: Failed to create updated topology files')
        print('Source ESP: {0:s}'.format(esp_filename))
        print('Source AC: {0:s}'.format(ac_filename))
        print('Constraints: {0:s}'.format(constraint_filename))
        sys.exit(1)

    return


def update_topology_for_combining(average_charges, mol_info, atom_info,
                                  src_dir, output_dir, role='initial'):
    """
    Create updated topology for molecule, constraining matched atoms to the
    average charge between the initial and final molecules.

    Args:
        average_charges (dict): Average charges (initial and final molecules)
                               for matched atoms
        mol_info (MolInfo): Topology, structure and filename information
                            for the molecule
        atom_info (dict): Atomic charge and naming information by index
        src_dir (str): Path to directory containing files used in original
                       topology creation
        output_dir (str): Destination path for output
        role (str): Select either 'initial' or 'final' molecule
                    (determines file and residue naming)

    Returns:
        str: Path to directory containing output topology files
        str: Path to output Amber library file

    """


    struct = mol_info.struct
    ac_filename = mol_info.ac_filename

    top_output_dir = os.path.join(output_dir, role+'_top')
    os.makedirs(top_output_dir)

    start_dir = os.getcwd()

    os.chdir(top_output_dir)

    # Write constraint file for use in updating topology
    constraint_filename = os.path.join(top_output_dir, 'chg_constraint')
    write_charge_constraint_file(average_charges, struct, atom_info,
                                 constraint_filename)

    if role == 'initial':
        mol_name = 'ini'
    else:
        mol_name = 'fin'

    # Create updated topology using charge constraints for matched atoms
    create_updated_prep_frcmod(src_dir, ac_filename, constraint_filename,
                               mol_name)

    lib_filename = prep_to_lib(mol_name)
    lib_filename = os.path.abspath(lib_filename)

    os.chdir(start_dir)

    return top_output_dir, lib_filename


def prep_to_lib(mol_name):
    """
    Convert an Amber prep file into an Amber library file (which is easier to
    manipulate).

    Args:
        mol_name (str): Three letter name used for molecules residue name

    Returns:
        str: Filename of output library file

    """

    prep_to_lib_template = '''
source leaprc.gaff
loadamberparams {0:s}.frcmod
loadamberprep {0:s}.prep
saveamberparm {1:s} {0:s}.top {0:s}.crd
saveoff {1:s} {0:s}.lib
quit
'''

    # Write tleap command to input file
    tleap_script = prep_to_lib_template.format(mol_name, mol_name.upper())

    leap_in_filename = mol_name+'.leapin'

    leap_in_file = open(leap_in_filename, 'w')
    print(tleap_script, file=leap_in_file)
    leap_in_file.close()

    # run tleap using input file
    returncode = subprocess.call(['tleap', '-s', '-f', leap_in_filename])

    if returncode:
        cur_dir = os.getcwd()
        print('ERROR: tleap failure in: {0:s}'.format(cur_dir))
        sys.exit(1)

    return '{0:s}.lib'.format(mol_name)


def create_combined_structure(initial_struct, final_struct,
                              common_names, resname='LIG'):
    """
    Create a sasmol.SasMol describing the hybrid ligand - i.e. copy the
    appearing atoms from the final molecule and add to those from the initial
    molecule.

    Args:
        initial_struct (sasmol.SasMol): Atomic information for the initial
                                        molecule
        final_struct (sasmol.SasMol): Atomic information for the final molecule
        common_names (list): Names of atoms common to both molecules
        resname (str): Name to be given to the combined molecule

    Returns:
        sasmol.SasMol:

    """

    frame = 0

    # Create a sasmol style selection string for the matched/common atoms
    formatted_atom_names = seg_list = ['"{0:s}"'.format(x) for x in common_names]
    filter_txt = 'name[i] not in [{0:s}]'.format(','.join(formatted_atom_names))

    # Select matched atoms in the final molecule and copy to new SasMol object
    err, mask = final_struct.get_subset_mask(filter_txt)

    if err:
        raise Exception(err)

    tmp_struct = sasmol.SasMol(0)
    final_struct.copy_molecule_using_mask(tmp_struct, mask, frame)

    # Combine the full ligand with disappearing atoms with just the appearing
    # atoms from the second ligand in a new SasMol object
    combined_structure = sasmol.SasMol(0)
    error = combined_structure.merge_two_molecules(initial_struct, tmp_struct)

    if error:
        raise Exception(error)

    # Tidy up the residue naming, etc. of the combined ligand
    natoms_combined = combined_structure.natoms()
    combined_structure.setResname([resname]*natoms_combined)
    combined_structure.setChain(['X']*natoms_combined)
    combined_structure.setSegname(['X']* natoms_combined)

    return combined_structure


def create_combined_amber_lib_file(initial_lib_filename, final_lib_filename,
                                   output_filename, common_atom_names,
                                   resname='LIG'):
    """
    Combine the Amber library files for the initial and final molecule to
    create one for the hybrid ligand.

    Args:
        initial_lib_filename (str): Path to amber library file for initial
                                    molecule
        final_lib_filename (str): Path to amber library file for final molecule
        output_filename (str): Path for output hybrid library file
        common_atom_names (list): Atom names common to both molecules
        resname (str): Name to use for hybrid residue

    Returns:

    """

    # Get atom and connectivity information from library files for both
    # molecules
    initial_lib_info = parse_amber_lib_file(initial_lib_filename)
    final_lib_info = parse_amber_lib_file(final_lib_filename)

    # Library indices have no correlation with those in PDB
    # Need to create new mapping to names
    initial_lib_match_idx_map = {}

    for idx in range(len(initial_lib_info['atom'])):

        atom = initial_lib_info['atom'][idx]

        if atom['name'] in common_atom_names:
            initial_lib_match_idx_map[atom['name']] = atom['index']


    n_atoms_original = len(initial_lib_info['atom'])

    # Edit library file information for the final molecule so that appearing
    # indices follow on from initial library file
    # Connections to matched section edited to point to equivalent atom in
    # initial molecule
    edited_final_lib_info = edit_final_lib_info(final_lib_info,
                                                n_atoms_original,
                                                initial_lib_match_idx_map)

    # Add appearing atoms from the final molecule to the library information
    initial_lib_info['connectivity'] += edited_final_lib_info['connectivity']
    initial_lib_info['connectivity'].sort(key=lambda x: x['ndx1'])

    initial_lib_info['atom'] += edited_final_lib_info['atom']

    initial_lib_info['positions'] += edited_final_lib_info['positions']

    write_amber_lib_from_info(initial_lib_info, output_filename, resname)

    return


def test_hybrid_library(output_dir):
    """
    Check that the output library can be used to create a valid amber topology.
    Also convert to prepc format anc heck with parmchk

    Args:
        output_dir (str): Path to directory holding output files

    Returns:

    """

    start_dir = os.getcwd()

    os.chdir(output_dir)

    test_system_build_template = '''
source leaprc.gaff
frcmod = loadamberparams hybrid.frcmod
loadoff hybrid.lib
hybrid = loadpdb hybrid.pdb
saveamberprep hybrid test.prepc
saveamberparm hybrid test.top test.crd
savepdb hybrid test.pdb
quit

'''

    leap_in_filename = 'test.leapin'

    leap_in_file = open(leap_in_filename, 'w')
    print(test_system_build_template, file=leap_in_file)
    leap_in_file.close()

    returncode = subprocess.call(['tleap', '-s', '-f', leap_in_filename])

    if returncode:
        print('ERROR: Test of the hybrid topology failed')
        sys.exit(1)

    returncode = subprocess.call(['parmchk',
                                  '-i', 'test.prepc',
                                  '-f', 'prepc', '-o', 'test.frcmod'])

    if returncode:
        print('ERROR: Unable to run parmchk on test.prepc')
        sys.exit(1)

    os.chdir(start_dir)

    return


def prepare_output_dir(output_path, delete_old):
    """
    Create output directory (and check it does not already exist)

    Args:
        output_path (str): Path to directory to store output
        delete_old (bool): Should existing directory be deleted then recreated

    Returns:

    """

    if os.path.exists(output_path):

        isdir = os.path.isdir(output_path)

        if isdir:

            print("Output directory already exists: {0:s}".format(output_path))

        else:

            print("Path for output directory is an existing file: {0:s}".format(output_path))

        if delete_old:

            if isdir:

                shutil.rmtree(output_path)

            else:

                os.remove(output_path)

        else:

            sys.exit(1)

    os.makedirs(output_path)
    os.makedirs(os.path.join(output_path, 'tmp'))

    return


def output_non_common_atom_names(structure, common_atom_names, output_filename):
    """
    Write the names of the atoms not in the common atoms list to a file.

    Args:
        structure (sasmol.SasMol): Atomic information about molecule
        common_atom_names (list): Name sof atoms in the common region
        output_filename (str): Path to output file

    Returns:

    """

    names = structure.name()

    uncommon_names = [x for x in names if x not in common_atom_names]

    out_txt = ' '.join(uncommon_names)

    out_file = open(output_filename, 'w')
    print(out_txt, file=out_file)
    out_file.close()

    return


def output_ties_input_files(structure, common_atom_names,
                            output_dir, role='initial'):
    """
    Ouput the files needed by BAC Builder to create TIES simulation input. The
    information output is different for initial and final molecules for
    historical reasons.

    Args:
        structure (sasmol.SasMol): Atomic information about molecule
        common_atom_names (list): Name sof atoms in the common region:
        output_dir (str): Path into which output files are saved
        role (str): Are we looking at the initial or final molecule?

    Returns:

    """

    if role == 'initial':
        output_filename = os.path.join(output_dir, 'namelist-lig1.dat')
    else:
        output_filename = os.path.join(output_dir, 'names-list.dat')

    names = structure.name()

    uncommon_names = [x for x in names if x not in common_atom_names]

    out_file = open(output_filename, 'w')

    for atom_name in uncommon_names:

        if role == 'initial':
            out_txt = atom_name
        else:
            # Placeholder used instead of original name
            out_txt = 'X\t{0:s}'.format(atom_name)

        print(out_txt, file=out_file)

    out_file.close()

    return


def parse_commandline_options():
    """
    Parse the commandline options and set some defaults

    Args:

    Returns:
        Namespace: Contains values parsed for each argument
    """

    parser = argparse.ArgumentParser(
        description='Combine two AMBER topologies (in prep format) for use in TIES')

    parser.add_argument(
        '-i', '--initial',
        required=True,
        help='Path to directory containing files describing the initial ligand')

    parser.add_argument(
        '-ip', '--ipdb',
        required=True,
        help='Path to PDB describing the initial ligand')

    parser.add_argument(
        '-f', '--final',
        required=True,
        help='Path to directory containing files describing the final ligand')

    parser.add_argument(
        '-fp', '--fpdb',
        required=True,
        help='Path to PDB describing the final ligand')

    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Path to directory into which to save output')

    parser.add_argument(
        '-t', '--task',
        choices=['compare', 'all'],
        default='all',
        help='Task to be performed')

    parser.add_argument(
        '-q', '--qtol',
        default=0.1,
        type=float,
        help='Charge difference tolerance for acceptable match substructure')

    parser.add_argument(
        '-a', '--atol',
        default=0.1,
        type=float,
        help='Charge difference tolerance for individual atoms in acceptable match substructure')

    parser.add_argument(
        '-m', '--manual',
        action='store_true',
        help='Flag to manually select the desired choice of common element')

    parser.add_argument(
        '-r', '--resname',
        default='LIG',
        help='Three letter residue name for hybrid molecule')

    parser.add_argument(
        '-n', '--noalign',
        action='store_true',
        help='Flag to prevent alignment of the ligands via the matched region (uses input coordinates)')

    parser.add_argument(
        '-d', '--delete_old',
        action='store_true',
        help='Flag to select deletion of directory to be created by the script before execution')

    parser.add_argument(
        '--amino',
        action='store_true',
        help='Flag to show that it should be assumed that the input is amino acid-like')

    args = parser.parse_args()

    return args


def main():
    """
    Run the comparison of two input ligands and create a hybrid PDB and library
    files for use in TIES simulations

    Returns:

    """

    args = parse_commandline_options()

    initial_dir = os.path.abspath(args.initial)
    initial_pdb = os.path.abspath(args.ipdb)
    final_dir = os.path.abspath(args.final)
    final_pdb = os.path.abspath(args.fpdb)
    output_dir = os.path.abspath(args.output)
    manual = args.manual
    q_tol = args.qtol
    q_atom_tol = args.atol
    align = not args.noalign
    output_resname = args.resname
    delete_old = args.delete_old
    amino = args.amino

    prepare_output_dir(args.output, delete_old)

    # Obtain structure and charge information on the two molecules
    # and get a selected common matched region
    # Submatch given and indices from the initial molecule
    (initial_mol_info,
     initial_atom_info,
     final_mol_info,
     final_atom_info,
     matched_idx_map,
     selected_submatch) = compare_ligands(initial_dir, initial_pdb,
                                          final_dir, final_pdb,
                                          output_dir, manual,
                                          tolerance=q_tol,
                                          atom_tolerance=q_atom_tol,
                                          amino=amino)

    if args.task == 'all':

        # Rename submatch atoms in the final molecule
        # This is propogated to updated PDB, prep and ac files
        update_final_description(final_mol_info, selected_submatch,
                                 initial_atom_info, matched_idx_map,
                                 output_dir)

        common_atom_names = [initial_atom_info[x].name for x in selected_submatch]

        if align:
            # Align final structure on initial
            align_molecules_using_common(final_mol_info.struct,
                                         initial_mol_info.struct,
                                         common_atom_names)

        # Calculate average charges for common atoms
        average_charges = calculate_average_charges(initial_atom_info,
                                                    final_atom_info,
                                                    selected_submatch,
                                                    matched_idx_map)

        # Prepare updated prep files using charge constraints from averages
        # Convert to lib files for easier processing

        (initial_top_dir,
         initial_libname) = update_topology_for_combining(average_charges,
                                                          initial_mol_info,
                                                          initial_atom_info,
                                                          initial_dir,
                                                          output_dir)

        (final_top_dir,
         final_libname) = update_topology_for_combining(average_charges,
                                                        final_mol_info,
                                                        final_atom_info,
                                                        final_dir,
                                                        output_dir,
                                                        role='final')

        # Create combined structure
        combined_structure = create_combined_structure(initial_mol_info.struct,
                                                       final_mol_info.struct,
                                                       common_atom_names,
                                                       resname=output_resname)

        hybrid_pdbname = os.path.join(output_dir, 'hybrid.pdb')
        combined_structure.write_pdb(hybrid_pdbname, 0, 'w')

        # Create combined lib
        hybrid_libname = os.path.join(output_dir, 'hybrid.lib')

        create_combined_amber_lib_file(initial_libname, final_libname,
                                       hybrid_libname, common_atom_names,
                                       resname=output_resname)

        # Create combined frcmod
        hybrid_frcmodname = os.path.join(output_dir, 'hybrid.frcmod')
        create_merged_frcmod(initial_mol_info.frcmod_filename,
                             final_mol_info.frcmod_filename,
                             hybrid_frcmodname)

        # Test hybrid ligand
        test_hybrid_library(output_dir)

        # Generate disappearing/appearing atom lists
        disappearing_atoms_filename = os.path.join(output_dir,
                                                   'disappearing_atoms.txt')

        output_non_common_atom_names(initial_mol_info.struct,
                                     common_atom_names,
                                     disappearing_atoms_filename)

        output_ties_input_files(initial_mol_info.struct,
                                common_atom_names,
                                output_dir,
                                role='initial')

        appearing_atoms_filename = os.path.join(output_dir,
                                                'appearing_atoms.txt')

        output_non_common_atom_names(final_mol_info.struct,
                                     common_atom_names,
                                     appearing_atoms_filename)

        output_ties_input_files(final_mol_info.struct,
                                common_atom_names,
                                output_dir,
                                role='final')

    return


if __name__ == "__main__":
    main()
