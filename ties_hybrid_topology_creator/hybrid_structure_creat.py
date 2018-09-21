import os
import math
import numpy as np
from sasmol.sasmol import SasMol
from rdkit import Chem

from bac_ties.bac_atom_data import *
from prepare_structure import prepare_param_for_matching


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
    tmp_struct = SasMol(0)
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

    align_coor, align_com = calc_coor_com(align_struct,common_atom_names)
    target_coor, target_com = calc_coor_com(target_struct, common_atom_names)
    msd = np.sum(np.square(np.subtract(align_coor, target_coor)))
    rmsd1 = math.sqrt(msd)

    align_struct.align(frame, align_coor, align_com, target_coor, target_com)

    align_coor, align_com = calc_coor_com(align_struct,common_atom_names)
    msd = np.sum(np.square(np.subtract(align_coor, target_coor)))
    rmsd2 = math.sqrt(msd)
    print ("NOTE: RMSDs of the selected MCS before and after alignment are %f	%f.\nIf the values too large, consider starting with better agreeing initial structures!!" % (rmsd1, rmsd2))

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
	float: Amount by which the overall charge of the common region 
		differs after averaging for the initial molecule
	float: Amount by which the overall charge of the common region 
		differs after averaging for the final molecule

    """

    average_charges = {}
    q_tot_initial = 0
    q_tot_final = 0

    for initial_idx in submatch:
        atom_name = initial_atom_info[initial_idx].name
        final_idx = matched_idx_map[initial_idx]

        charge1 = initial_atom_info[initial_idx].charge
	q_tot_initial += charge1
        charge2 = final_atom_info[final_idx].charge
	q_tot_final += charge2

        average_charges[atom_name] = (charge1 + charge2) / 2.0

    q_tot_updated = (q_tot_initial + q_tot_final) / 2.0
    del_q_initial = q_tot_updated - q_tot_initial    
    del_q_final = q_tot_updated - q_tot_final    

    return average_charges, del_q_initial, del_q_final


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
    formatted_atom_names = seg_list = [
        '"{0:s}"'.format(x) for x in common_names]
    filter_txt = 'name[i] not in [{0:s}]'.format(
        ','.join(formatted_atom_names))

    # Select matched atoms in the final molecule and copy to new SasMol object
    err, mask = final_struct.get_subset_mask(filter_txt)

    if err:
        raise Exception(err)

    tmp_struct = SasMol(0)
    final_struct.copy_molecule_using_mask(tmp_struct, mask, frame)

    # Combine the full ligand with disappearing atoms with just the appearing
    # atoms from the second ligand in a new SasMol object
    combined_structure = SasMol(0)
    error = combined_structure.merge_two_molecules(initial_struct, tmp_struct)

    if error:
        raise Exception(error)

    # Tidy up the residue naming, etc. of the combined ligand
    natoms_combined = combined_structure.natoms()
    combined_structure.setResname([resname] * natoms_combined)
    combined_structure.setResid([1] * natoms_combined)
    combined_structure.setChain(['X'] * natoms_combined)
    combined_structure.setSegname(['X'] * natoms_combined)

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
                             initial_atom_info, idx_map, output_dir, atom_type):
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
	atom_type (str): Amber atom type

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
                               renamed_pdb, renamed_ac, renamed_prep, atom_type)

    os.chdir(cwd)

    final_mol_info.ac_filename = renamed_ac
    final_mol_info.prep_filename = renamed_prep
    final_mol_info.pdb_filename = renamed_pdb

    return
