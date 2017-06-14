### Workflow documentation:

1.  Read in global parameters from user:
    Inital and Final PDB and paths
    Output directory
    Task, charge tolerance, manual or auto sss
    Resname (LIG), DON'T align, delete old in dir
    Amino acid like?

2.  Prepare output directories. Deletes old stuff in there

3.  Compare ligands (detailed info at 3.X)
    Returns ???

4.  Update stuff in the final files
	4.1 Replace the name of atoms in the final structure based on
	index, with the name from the initial file.
	4.2 Covert ac and prep files of the final to
	match the new naming convention. Antechamber is used.

5.	Align (move) the final structure coordinates to match
	the inital structure coordinates based on the common substructure.
	A lot of SASMOL is used here. Uhh.

6. 	Calculate the average charges in the common atoms list. Nothing fancy.

7. 	This does a couple of thing written down in .X
	7.1 Basically just writes a file for antechamber to contrain charges.
	7.2 Antechamber is used again for some kind of transformations.
	7.3 Convert prep to lib. LEAP is used.

8. 	Paste the appearing part into the inital strucutre.



### Call tree:

1.  parse_commandline_options							-> parse.py
2.  prepare_output_dir									-> directory.py
3.  compare_ligands                                     -> compare.py
	1. prepare_mols 									-> prepare.py
		1. get_param_files								-> prepare.py
			1. get_user_selection						-> user_selection.py
		2. prepare_mol_for_matching						-> prepare.py
			1. prepare_structure_for_matching			-> prepare.py
				1. create_safe_atom_name_pdb			-> prepare.py
					1.  element_from_pdb_atom_name		-> prepare.py
			2. prepare_param_for_matching				-> prepare.py
		3. get_param_files								-> prepare.py
		4. prepare_mol_for_matching						-> prepare.py
	2. create_charge_idx_map                            -> bac_atom_data.py
	3. get_atom_info_rdkit                              -> bac_atom_data.py
	4. remove_ring_error                                -> compare.py
	5. check_atom_type_match                            -> compare.py
	6. get_submatches                                   -> compare.py
		1. get_bridge_atoms                             -> compare.py
		2. get_stop_ring_remove                         -> compare.py
			1. get_bridge_ring                          -> compare.py
		3. get_linked_to_remove                         -> compare.py
	7. get_charge_diffs_submatches                      -> compare.py
	8. get_user_submatch                                -> user_selection.py
		1. get_user_selection							-> user_selection.py
	9. output_submatches_file     						-> output.py
4.  update_final_description
    	1. rename_common_atoms_final
		2. prepare_param_for_matching
5.  align_molecules_using_common   						-> structure_align.py
		1. calc_coor_com		  						-> structure_align.py
		2. align(SASMOL)		   						-> SASMOL
6.  calculate_average_charges      						-> hybrid_create.py
7.  update_topology_for_combining       				-> topology.py
		1. write_charge_constraint_file 				-> topology.py
		2. create_updated_prep_frcmod   				-> topology.py
		3. prep_to_lib 									-> topology.py
8.  create_combined_structure	   						-> hybrid_create.py
9.  create_combined_amber_lib_file 						-> hybrid_create.py
10. create_merged_frcmod           						-> bac_amber_utils.py
11. test_hybrid_library            						-> output.py
12. output_non_common_atom_names   						-> output.py
13. output_ties_input_files        						-> output.py
