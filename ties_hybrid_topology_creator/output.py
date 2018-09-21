from __future__ import print_function

import os
import subprocess
import sys
import shutil
from stat import *

def output_MCS_file(initial_info, final_info,
                           matched_idx_map, output_dir):
    """
    Output a summary file detailing the possible submatches and the atomic
    charge differences between the initial and final molecules

    Args:
        initial_info (dict): Naming and connectivity information (AtomInfo)
                             by atom index (initial molecule)
        final_info (dict): Naming and connectivity information (AtomInfo)
                           by atom index (final molecule)
        matched_idx_map (dict): Map index of atom in initial to that in the
                                final molecule
        output_dir (str): Path in which to save output

    Returns:

    """

    out_filename = os.path.join(output_dir, 'MCS_details.txt')

    out_file = open(out_filename, 'w')

    print('Uses atom names from initial molecule\n', file=out_file)
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
    unmatched_idxs = [x for x in initial_info.keys() if x not in matched_idx_map.keys()]
    print('Name\tQ', file=out_file)
    for idx in unmatched_idxs:
        name = initial_info[idx].name
	charge = initial_info[idx].charge
        print('{0:s}\t{1:f}'.format(name, charge),
              file=out_file)

    print('\nUnmatched atoms from final:', file=out_file)
    unmatched_idxs = [x for x in final_info.keys(
    ) if x not in matched_idx_map.values()]
    print('Name\tQ', file=out_file)
    for idx in unmatched_idxs:
        name = final_info[idx].name
	charge = final_info[idx].charge
        print('{0:s}\t{1:f}'.format(name, charge),
              file=out_file)

    out_file.close()

    return

def output_non_common_atom_names(
        structure,
        common_atom_names,
        output_filename):
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


def test_hybrid_library(output_dir, atom_type):
    """
    Check that the output library can be used to create a valid amber topology.
    Also convert to prepc format anc heck with parmchk

    Args:
        output_dir (str): Path to directory holding output files
	atom_type (str): AMBER atom type

    Returns:

    """

    start_dir = os.getcwd()

    os.chdir(output_dir)

    test_system_build_template = '''
source leaprc.{0:s}
frcmod = loadamberparams hybrid.frcmod
loadoff hybrid.lib
hybrid = loadpdb hybrid.pdb
saveamberprep hybrid test.prepc
saveamberparm hybrid test.top test.crd
savepdb hybrid test.pdb
quit

'''

    test_system_build_script = test_system_build_template.format(atom_type)

    leap_in_filename = 'test.leapin'

    leap_in_file = open(leap_in_filename, 'w')
    print(test_system_build_script, file=leap_in_file)
    leap_in_file.close()

    output = subprocess.check_output(
        ['tleap', '-s', '-f', leap_in_filename], stderr=subprocess.STDOUT)

    missing_angles = []
    missing_dihedrals = []

    for line in output.splitlines():

        if "Could not find angle parameter:" in line:

            cols = line.split(':')
            angle = cols[-1].strip()
            if angle not in missing_angles:
                missing_angles.append(angle)

        elif "No torsion terms for" in line:
	    cols = line.split()
	    torsion = cols[-1].strip()
            if torsion not in missing_dihedrals:
                missing_dihedrals.append(torsion)

    if missing_angles or missing_dihedrals:

        print('WARNING: Adding default values for missing dihedral to frcmod')
#        print('WARNING: Okay unless there are atom type changes in match')

        old_frcmod = open('hybrid.frcmod')
        frcmod_lines = old_frcmod.readlines()
        old_frcmod.close()

        new_frcmod = open('hybrid.frcmod', 'w')

        for line in frcmod_lines:

            new_frcmod.write(line)

            if 'ANGLE' in line:
                for angle in missing_angles:
                    new_frcmod.write(
                        '{:<14}0     120.010   same as ca-ca-ha\n'.format(angle))

            if 'DIHE' in line:
                for angle in missing_dihedrals:
                    new_frcmod.write(
                        '{:<14}1    0.00       180.000           2.000      same as X -c2-ca-X\n'.format(angle))

	new_frcmod.close()

        returncode = subprocess.call(['tleap', '-s', '-f', leap_in_filename])

        if returncode:
            print('ERROR: Test of the hybrid topology failed')
            sys.exit(1)

    if not os.path.exists('test.top'):

        print('ERROR: Unable to create test.top, running parmchk')

	if atom_type == 'gaff':
	    at = 1
	else:
	    at = 2

        returncode = subprocess.call(['parmchk2',
                                      '-i', 'test.prepc',
                                      '-f', 'prepc', '-o', 'missing.frcmod', '-s', str(at)])

        if returncode:
            print('ERROR: Unable to run parmchk on test.prepc')
            sys.exit(1)

    if os.stat('test.top').st_size == 0:
        subprocess.call(['tleap', '-s', '-f', 'test.leapin'], stderr=subprocess.STDOUT)    

    os.chdir(start_dir)

    print('\nHybrid topology created correctly')

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

            print(
                "Path for output directory is an existing file: {0:s}".format(output_path))

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
