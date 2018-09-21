from __future__ import print_function

import os
import shutil
import subprocess
import sys

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

            charge_diff = abs(avg_charge - original_charge)

            if charge_diff > 0.5:
                print(
                    "Charge difference > 0.5 (0:f) for atom {1:s}".format(
                        charge_diff, atom_name))

            print('CHARGE {0:.6f} {1:d} {2:s}'.format(avg_charge,
                                                      idx + 1,
                                                      atom_name),
                  file=out_file)

    out_file.close()

    return

def create_updated_prep_frcmod(src_dir, ac_filename,
                               constraint_filename, mol_name, atom_type):
    """
    Use antechamber to produce prep and frcmod files for selected molecule,
    constraining selected atomic charges to average values (between initial
    and final molecules).

    Args:
        src_dir (str): Path in which to find ESP file
        ac_filename (str): Path to original ac file
        constraint_filename (str): Path to charge constraint file
        mol_name (str): Name to be given to output molecule
	atom_type (str): Amber atom type

    Returns:

    """

    source_esp_filename = os.path.join(src_dir, 'ANTECHAMBER.ESP')
    esp_filename = 'original_antrechamber.esp'
    shutil.copyfile(source_esp_filename, esp_filename)

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
                                   '-rn', mol_name.upper(),
				   '-at', atom_type])

    # Output missing force field parameters in frcmod format
    if atom_type == 'gaff':
	at = 1
    else:
	at = 2

    returncode += subprocess.call(['parmchk2', '-i', mol_name + '.prep',
                                   '-f', 'prepi',
                                   '-o', mol_name + '.frcmod', '-s', str(at)])

    if returncode:
        print('ERROR: Failed to create updated topology files')
        print('Source ESP: {0:s}'.format(esp_filename))
        print('Source AC: {0:s}'.format(ac_filename))
        print('Constraints: {0:s}'.format(constraint_filename))
        sys.exit(1)

    return


def create_updated_prep_frcmod_for_others(ac_filename,
                               constraint_filename, mol_name, atom_type):
    """
    Use antechamber to produce prep and frcmod files for selected molecule,
    constraining all atomic charges to pre-defined values.

    Args:
        ac_filename (str): Path to original ac file
        constraint_filename (str): Path to charge constraint file
        mol_name (str): Name to be given to output molecule
	atom_type (str): Amber atom type

    Returns:

    """

    returncode = 0

    # Get input files for two stage resp fitting (using constraints).
    # Here, we need not perform resp calculation, so this step only
    # helps us prepare the charge file to be read in while writing out prep file.
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
    # Remove redundant resp input files
    os.remove(mol_name + '.respin1')
    os.remove(mol_name + '.respin2')

    # Convert ac file into prep but reading in charges taken from the constraint file
    returncode += subprocess.call(['antechamber', '-i', ac_filename,
                                   '-fi', 'ac', '-c', 'rc',
                                   '-cf', qin_filename,
                                   '-o', mol_name + '.prep',
                                   '-fo', 'prepi',
                                   '-rn', mol_name.upper(),
				   '-at', atom_type])

    # Output missing force field parameters in frcmod format
    if atom_type == 'gaff':
	at = 1
    else:
	at = 2

    returncode += subprocess.call(['parmchk2', '-i', mol_name + '.prep',
                                   '-f', 'prepi',
                                   '-o', mol_name + '.frcmod', '-s', str(at)])

    if returncode:
        print('ERROR: Failed to create updated topology files')
        print('Source AC: {0:s}'.format(ac_filename))
        print('Constraints: {0:s}'.format(constraint_filename))
        sys.exit(1)

    return


def update_topology_for_combining(average_charges, mol_info, atom_info,
                                  src_dir, output_dir, del_q, charge, num_common_atoms, atom_type, role='initial'):
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
	del_q (float): Amount by which the overall charge of the matched
		    region differs from its original value.
	charge (str): Charge method taken as an argument (resp or others)
	num_common_atoms (int): Number of atoms in the matched region
	atom_type (str): Amber atom type
        role (str): Select either 'initial' or 'final' molecule
                    (determines file and residue naming)
    Returns:
        str: Path to directory containing output topology files
        str: Path to output Amber library file

    """

    struct = mol_info.struct
    ac_filename = mol_info.ac_filename

    top_output_dir = os.path.join(output_dir, role + '_top')
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

    # If the charge method is anything other than RESP, calculate the new charges
    # for the atoms in unmatched region by uniformly distributing del_q
    # over the unmatched region of the input molecule. Thereafter, append
    # the new charges to the constraint file such that charges of all atoms are constrained to
    # their new values before creating updated prep/frcmod files

    # Create updated topology using charge constraints for matched atoms (in the case of RESP charges)
    # or all atoms (in the case other charge methods)

    if charge == 'others':
	append_charge_constraint_file(constraint_filename, average_charges, struct, atom_info, del_q, num_common_atoms)
	create_updated_prep_frcmod_for_others(ac_filename, constraint_filename,
                               mol_name, atom_type)
    else:
	create_updated_prep_frcmod(src_dir, ac_filename, constraint_filename,
                               mol_name, atom_type)

    lib_filename = prep_to_lib(mol_name, atom_type)
    lib_filename = os.path.abspath(lib_filename)

    os.chdir(start_dir)

    return top_output_dir, lib_filename


def append_charge_constraint_file(constraint_filename, avg_charges, struct, atom_info, del_q, num_common_atoms):
    """
    Calculates the new charges for the atoms in unmatched region by uniformly distributing
    del_q over the unmatched region of the input molecule. 
    Thereafter, append the new charges to the constraint file such that charges of all atoms 
    are constrained to their new values before creating updated prep/frcmod files using 
    constrained RESP calculation.

    Args:
	constraint_filename (str): Path to the constraint filename.
        avg_charges (dict): Average charges (initial and final molecules)
                               for matched atoms
        struct (sasmol.SasMol: Atomic information for molecule
	atom_info (dict): Atom charge and naming info by index
	del_q (float): Amount by which the overall charge of the matched
		    region differs from its original value.
	num_common_atoms (int): Number of atoms in the matched region.

    Returns:

    """
    tot_num_atoms = struct.natoms()
    num_unmatched_atoms = tot_num_atoms - num_common_atoms
    del_q_per_atom = del_q/num_unmatched_atoms

    out_file = open(constraint_filename, 'a')

    for idx, atom_name in enumerate(struct.name()):

        if atom_name not in avg_charges:

            # idx use here is a bit dubious - but should work
            original_charge = atom_info[idx].charge
	    new_charge = original_charge - del_q_per_atom 

            charge_diff = abs(new_charge - original_charge)

            if charge_diff > 0.5:
                print(
                    "Charge difference > 0.5 (0:f) for atom {1:s}".format(
                        charge_diff, atom_name))

            print('CHARGE {0:.6f} {1:d} {2:s}'.format(new_charge,
                                                      idx + 1,
                                                      atom_name),
                  file=out_file)

    out_file.close()

    return


def prep_to_lib(mol_name, atom_type):
    """
    Convert an Amber prep file into an Amber library file (which is easier to
    manipulate).

    Args:
        mol_name (str): Three letter name used for molecules residue name
	atom_type (str): AMBER atom type.

    Returns:
        str: Filename of output library file

    """

    prep_to_lib_template = '''
source leaprc.{2:s}
loadamberparams {0:s}.frcmod
loadamberprep {0:s}.prep
saveamberparm {1:s} {0:s}.top {0:s}.crd
saveoff {1:s} {0:s}.lib
quit
'''
    # Write tleap command to input file
    tleap_script = prep_to_lib_template.format(mol_name, mol_name.upper(), atom_type)

    leap_in_filename = mol_name + '.leapin'

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

