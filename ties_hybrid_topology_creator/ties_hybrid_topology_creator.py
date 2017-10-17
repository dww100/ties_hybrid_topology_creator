#! /usr/bin/env python
"""
Build a hybrid Amber topology from single topologies for two ligands ready for
use in TI calculations.
Elements of the initial molecule 'disappear' whilst those of the final 'appear'
"""
from __future__ import print_function

import argparse

from compare_structure import compare_ligands
from hybrid_structure_creat import *
from topology import update_topology_for_combining
from output import *


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
        '-a',
        '--atol',
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
