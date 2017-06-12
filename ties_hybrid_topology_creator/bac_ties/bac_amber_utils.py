#! /usr/bin/env python
"""

"""
from __future__ import print_function

import sasmol.sasmol as sasmol

from bac_txt_utils import *


class AmberPrep():
    def __init__(self, filename=''):

        self.filename = filename
        self.atom_name_from_idx = {}
        self.atom_type_from_name = {}
        self.h_links = {}
        self.resname = ''

        if filename:
            self.read_prep(filename)

    def read_prep(self, filename):
        """
        Extract the residue name and atom names of molecule described by a
        prep file.

        Args:
            filename (str): Path to prep file

        Returns:
            NA

        """

        idx_name_map = self.atom_name_from_idx
        h_idx_link = {}

        with open(filename) as prep_file:

            # Ignore first 4 header lines
            for _ in range(4):
                prep_file.readline()

            # Readline that describes residue contents
            line = prep_file.readline()
            self.resname = line.split()[0]

            # Skip remainder of header + dummy atom lines
            for _ in range(5):
                prep_file.readline()

            for line in prep_file:

                # Stop when blank line is reached (end of section)
                if not line.strip():
                    break

                cols = line.split()

                atom_idx = cols[0]

                # Amber likes to change CL -> Cl and BR -> Br
                atom_name = cols[1].upper()

                atom_type = cols[2]

                connection_idx = cols[4]

                idx_name_map[atom_idx] = atom_name

                self.atom_type_from_name[atom_name] = atom_type

                if atom_name.strip('0123456789') == 'H':
                    h_idx_link[atom_idx] = connection_idx

        for h_idx, connection_idx in h_idx_link.iteritems():
            h_name = idx_name_map[h_idx]
            connection_name = idx_name_map[connection_idx]
            self.h_links[h_name] = connection_name

        return


def get_resp_charge_per_atom_from_ac(ac_filename, structure):
    """
    Read in atomic charges from antechamber resp file and return dictionary
    mapping them to the atom names in the structure.

    Note: There are atom names in ac files but these may not agree with those
          read from the PDB

    Args:
        ac_filename (str): Path to the ac format file containing resp charges
        structure (sasmol.SasMol): Coordinate and naming information from PDB

    Returns:
        dict: Mapping of atom name to charge
    """

    charges = []

    structure_atom_names = structure.name()

    idx = 0

    with open(ac_filename) as ac_file:

        for line in ac_file:

            if line.startswith('ATOM'):

                cols = line.split()
                charges.append(float(cols[8]))

                atom_index = cols[1].strip()
                atom_name = cols[2]

                if remove_digits(atom_name).upper() != remove_digits(
                        structure_atom_names[idx]):
                    remove_digits(atom_name)

                    err_text = "Element mismatch between {0:s} and input PDB for atom {1:s}: {2:s} vs {3:s}".format(
                        ac_filename, atom_index, atom_name,
                        structure_atom_names[idx])
                    raise Exception(err_text)

                idx += 1

    if len(structure_atom_names) == len(charges):

        atom_name_charge_map = dict(zip(structure_atom_names, charges))

    else:

        raise Exception(
            "No. charges in ac (resp charge) file does not match no. atoms in PDB")

    return atom_name_charge_map


def read_pdb_prep_pair(prep_filename, pdb_filename):
    """
    Load structure, topology and naming information from PDB and atomic
    charges from antechamber ac file.

    Args:
        prep_filename (str): Path to the library (prep) file
        pdb_filename (str): Path to the PDB containing atomic coordinates
                            and names

    Returns:
        amber_utils.AmberPrep: Information derived from prep file
        sasmol.SasMol: Information read from PDB file

    """

    prep = AmberPrep(prep_filename)

    structure = sasmol.SasMol(0)
    structure.read_pdb(pdb_filename)
    structure.setFilename(pdb_filename)

    if len(structure.resnames()) != 1:
        raise Exception("More than one residue in PDB: {0:s}".format(pdb_filename))

    # Sometimes building ligands files end up with inconsistent atom name
    # capitilization - we will upper case everything

    upper_names = []
    structure_names = structure.name()

    for idx in range(structure.natoms()):
        upper_names.append(structure_names[idx].upper())

    structure.setName(upper_names)

    invalid = check_prep_structure_consistency(prep, structure)

    if invalid:
        raise Exception(invalid)

    return prep, structure


def check_prep_structure_consistency(prep, structure):
    """
    Check if the library (prep) file and structure read from PDB contain
    information on the same atoms.

    Args:
        prep (amber_utils.AmberPrep): Amber prep file information
        structure (sasmol.SasMol): PDB information

    Returns:
        bool: Do the prep and PDB contain the same atoms
    """

    err = ''

    file_text = " {0:s} {1:s}".format(prep.filename, structure.filename())

    if structure.resnames()[0] != prep.resname:
        err = "ERROR: Resname mismatch between PDB and prep file -" + file_text

    structure_names = structure.name()

    prep_names = prep.atom_name_from_idx.values()
    prep_h_names = prep.h_links.keys()

    if len(prep.atom_name_from_idx) != structure.natoms():
        err = "ERROR: Atom number mismatch -" + file_text

    # Hydrogens can have odd issues with name changes in prep so don't compare
    elif not (set(prep_names) - set(prep_h_names)).issubset(set(structure_names)):

        err = "ERROR: Heavy atom mismatch -" + file_text

    return err


def parse_amber_lib_file(filename):
    lib_info = {}

    sections = {'!!index array': 'ignore',
                'unit.atoms table': 'atom',
                'unit.atomspertinfo table': 'ignore',
                'unit.boundbox': 'ignore',
                'unit.childsequence': 'ignore',
                'unit.connect array': 'ignore',
                'unit.connectivity table': 'connectivity',
                'unit.hierarchy table': 'ignore',
                'unit.name single': 'ignore',
                'unit.positions table': 'positions',
                'unit.residueconnect table': 'ignore',
                'unit.residues table': 'ignore',
                'unit.residuesPdbSequenceNumber': 'ignore',
                'unit.solventcap': 'ignore',
                'unit.velocities': 'ignore'}

    section_heading_patterns = sections.keys()

    section = 'ignore'
    with open(filename) as f:

        for line in f:

            if line.startswith('!'):
                for pattern, process_choice in sections.iteritems():
                    if pattern in line:
                        section = process_choice
                        if section != 'ignore':
                            lib_info[section] = []

            elif section != 'ignore':

                if section == 'atom':

                    cols = line.split()
                    info = {}
                    info['type_txt'] = ' '.join(cols[1:5])
                    info['property_txt'] = ' '.join(cols[6:])
                    info['index'] = int(cols[5])
                    info['name'] = cols[0].strip('"').upper()
                    lib_info['atom'].append(info)

                elif section == 'connectivity':
                    cols = line.split()
                    info = {'ndx1': int(cols[0]),
                            'ndx2': int(cols[1])}
                    lib_info['connectivity'].append(info)

                elif section == 'positions':
                    lib_info['positions'].append(line)

    return lib_info


def edit_final_lib_info(lib_info, last_dis_ndx, matched_name_to_original_idxs):
    appearing_count = 0
    new_index_map = {}

    edited_info = {'atom': [],
                   'connectivity': [],
                   'positions': []}

    lib_idx_name_map = {}

    for i in range(len(lib_info['atom'])):

        atom = lib_info['atom'][i]

        if atom['name'] not in matched_name_to_original_idxs:

            # Indices of appearing part start after disappearing atoms
            appearing_count += 1
            new_index = last_dis_ndx + appearing_count
            new_index_map[atom['index']] = new_index

            info = atom
            info['index'] = new_index
            edited_info['atom'].append(info)

            # positions ordered as the atom information
            edited_info['positions'].append(lib_info['positions'][i])

        else:

            lib_idx_name_map[atom['index']] = atom['name']

    for bond in lib_info['connectivity']:

        if (bond['ndx1'] in new_index_map) and (bond['ndx2'] in new_index_map):

            new_bond = {'ndx1': new_index_map[bond['ndx1']],
                        'ndx2': new_index_map[bond['ndx2']]}
            edited_info['connectivity'].append(new_bond)

        elif (bond['ndx1'] in lib_idx_name_map) and (bond['ndx2'] not in lib_idx_name_map):

            name = lib_idx_name_map[bond['ndx1']]

            new_bond = {'ndx1': matched_name_to_original_idxs[name],
                        'ndx2': new_index_map[bond['ndx2']]}

            edited_info['connectivity'].append(new_bond)

        elif (bond['ndx1'] not in lib_idx_name_map) and (bond['ndx2'] in lib_idx_name_map):

            name = lib_idx_name_map[bond['ndx2']]

            new_bond = {'ndx1': matched_name_to_original_idxs[name],
                        'ndx2': new_index_map[bond['ndx1']]}

            edited_info['connectivity'].append(new_bond)

    return edited_info


def write_amber_lib_from_info(lib_info, output_filename, ligname):
    n_atoms = len(lib_info['atom'])

    out_file = open(output_filename, 'w')

    out_file.write('!!index array str\n')
    out_file.write(' "{0:s}"\n'.format(ligname))

    out_file.write('!entry.{0:s}.unit.atoms table  str name  str type  '
                   'int typex  int resx  int flags  int seq  int elmnt  '
                   'dbl chg\n'.format(ligname))
    for atom in lib_info['atom']:
        line_txt = ' "{0:s}" {1:s} {2:d} {3:s}\n'.format(atom['name'],
                                                         atom['type_txt'],
                                                         atom['index'],
                                                         atom['property_txt'])
        out_file.write(line_txt)

    out_file.write('!entry.{0:s}.unit.atomspertinfo table  str pname  '
                   'str ptype  int ptypex  int pelmnt  dbl pchg\n'.format(ligname))
    for atom in lib_info['atom']:
        atom_type = atom['type_txt'].split()[0]
        line_txt = ' "{0:s}" {1:s} 0 -1 0.0\n'.format(atom['name'], atom_type)
        out_file.write(line_txt)

    std_lines = '''!entry.{0:s}.unit.boundbox array dbl
 -1.000000
 0.0
 0.0
 0.0
 0.0
!entry.{0:s}.unit.childsequence single int
 2
!entry.{0:s}.unit.connect array int
 1
 {1:d}
'''.format(ligname, n_atoms)

    out_file.write(std_lines)

    out_file.write('!entry.{0:s}.unit.connectivity table  int atom1x  int atom2x  int flags\n'.format(ligname))
    for bond in lib_info['connectivity']:
        line_txt = ' {0:d} {1:d} 1\n'.format(bond['ndx1'],
                                             bond['ndx2'])
        out_file.write(line_txt)

    out_file.write('!entry.{0:s}.unit.hierarchy table  str abovetype  '
                   'int abovex  str belowtype  int belowx\n'.format(ligname))
    out_file.write(' "U" 0 "R" 1\n')
    for i in range(1, n_atoms + 1):
        out_file.write(' "R" 1 "A" {0:d}\n'.format(i))

    out_file.write('!entry.{0:s}.unit.name single str\n'.format(ligname))
    out_file.write(' "{0:s}"\n'.format(ligname))

    out_file.write('!entry.{0:s}.unit.positions table  dbl x  dbl y  dbl z\n'.format(ligname))
    for line in lib_info['positions']:
        out_file.write(line)

    std_lines = '''!entry.{0:s}.unit.residueconnect table  int c1x  int c2x  int c3x  int c4x  int c5x  int c6x
 1 {1:d} 0 0 0 0
!entry.{0:s}.unit.residues table  str name  int seq  int childseq  int startatomx  str restype  int imagingx
 "{0:s}" 1 {2:d} 1 "?" 0
!entry.{0:s}.unit.residuesPdbSequenceNumber array int
 0
!entry.{0:s}.unit.solventcap array dbl
 -1.000000
 0.0
 0.0
 0.0
 0.0
'''.format(ligname, n_atoms, n_atoms + 1)

    out_file.write(std_lines)

    out_file.write('!entry.{0:s}.unit.velocities table  dbl x  dbl y  dbl z\n'.format(ligname))
    for i in range(n_atoms):
        out_file.write(' 0.0 0.0 0.0\n')

    out_file.close()

    return


def parse_frcmod_sections(filename):
    frcmod_info = {}
    section = 'REMARK'

    with open(filename) as f:

        for line in f:

            start_line = line[0:9].strip()

            if start_line in ['MASS', 'BOND', 'IMPROPER',
                              'NONBON', 'ANGLE', 'DIHE']:
                section = start_line
                frcmod_info[section] = []

            elif line.strip() and section != 'REMARK':

                frcmod_info[section].append(line)

    return frcmod_info


def create_merged_frcmod(filename1, filename2, output_filename):
    frcmod_info1 = parse_frcmod_sections(filename1)
    frcmod_info2 = parse_frcmod_sections(filename2)

    output_file = open(output_filename, 'w')

    output_file.write('merged frcmod\n')

    for section in ['MASS', 'BOND', 'ANGLE',
                    'DIHE', 'IMPROPER', 'NONBON']:
        section_lines = set(frcmod_info1[section] + frcmod_info2[section])
        output_file.write('{0:s}\n'.format(section))
        for line in section_lines:
            output_file.write('{0:s}'.format(line))
        output_file.write('\n')

    output_file.write('\n')
    output_file.write('\n')

    output_file.close()

    return
