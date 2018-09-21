import itertools

from rdkit.Chem import rdFMCS

from prepare_structure import *
from bac_ties.bac_atom_data import *
from user_selection import *
from output import *


def compare_ligands(initial_dir, initial_pdb, final_dir, final_pdb,
                    output_dir, atom_type, tolerance=0.1, atom_tolerance=0.1,
                    Hatom_tolerance=0.1, amino=False):
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
	atom_type (str): Amber atom type

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
                                                    output_dir, atom_type, amino=amino)

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
    mcs_result = rdFMCS.FindMCS(rdkit_mols, completeRingsOnly=True, atomCompare=rdFMCS.AtomCompare.CompareElements, bondCompare=rdFMCS.BondCompare.CompareOrderExact, matchValences=True, ringMatchesRingOnly=True, timeout=60, maximizeBonds=False)

    if mcs_result:

        shared = Chem.MolFromSmarts(mcs_result.smartsString)

        initial_match_idxs = initial_mol_info.mol.GetSubstructMatch(shared)
        final_match_idxs = final_mol_info.mol.GetSubstructMatch(shared)

        # Map to get equivalent atom index in final molecule to the one in the
        # initial model (matched region)
        matched_idx_map = dict(zip(initial_match_idxs, final_match_idxs))

	
	# Delete all the H-atoms and heavy atoms (along with all the attached Hs) from the MCS
	# with atomic charges differing more than the atom_tolerance 
	initial_match_idxs = delete_atoms_with_large_qdiff(initial_match_idxs, matched_idx_map,
                                                         initial_charge_map, final_charge_map, 
							 initial_atom_info, atom_tolerance, Hatom_tolerance)

        # Update final matched idx list and initial to final map
        final_match_idxs = [matched_idx_map[x] for x in initial_match_idxs]
        matched_idx_map = dict(zip(initial_match_idxs, final_match_idxs))


	# Delete all Hs and heavy atoms (along with all attached Hs) with differing atom types from the MCS
	initial_match_idxs = correct_atom_type_mismatch(initial_match_idxs, matched_idx_map,
							initial_mol_info, final_mol_info, initial_atom_info)


	# Check if the the overall charge difference of the MCS between the initial and final
	# molecules is within the tolerance. If not, delete from MCS the atom or the group of atoms
	# (in case of heavy atom) with the largest q_diff. Repeat the process until the overall Q_diff is within the tolerance.
       	initial_match_idxs = overall_qdiff_correction(initial_match_idxs, matched_idx_map, initial_charge_map, 
							final_charge_map, initial_atom_info, tolerance)

	
	# Write out the initial molecule's structure in pdb format; 
	# Delete the atoms not in the matched region (to be used to find disjoint fragments in the MCS)
	tmp_pdbname = os.path.join(output_dir, 'temp.pdb')
	initial_mol_info.struct.write_pdb(tmp_pdbname, 0, 'w')

	tmp_pdbout = os.path.join(output_dir, 'temp_out.pdb')
	with open(tmp_pdbname) as infile, open(tmp_pdbout, 'w') as outfile:
	    	for line in infile:
	        	if any(' ' + initial_atom_info[idx].name + ' ' in line for idx in initial_match_idxs) or ('END' in line):
		            outfile.write(line)

	infile.close()
	outfile.close()

	# Check if there are any disjoint fragments in the MCS;
	# If yes, then take the largest of them as the final MCS
	matched_mol = find_disjoint_fragments(tmp_pdbout)
	
	if matched_mol != None:
		new_mol = Chem.MolFromSmarts(matched_mol)
		new_mol_idxs = initial_mol_info.mol.GetSubstructMatch(new_mol)
#		new_initial_match_idxs = []
#		for idx in new_mol_idxs:
#			if idx in initial_match_idxs:
#				new_initial_match_idxs.append(idx)
#				for bound in initial_atom_info[idx].bound:
#		                      if (initial_atom_info[bound].name.strip('0123456789') == 'H') and (bound in initial_match_idxs):
#	        	                new_initial_match_idxs.append(bound)
#				
#		initial_match_idxs = new_initial_match_idxs
		initial_match_idxs = new_mol_idxs
	
        # Update final matched idx list and initial to final map

        final_match_idxs = [matched_idx_map[x] for x in initial_match_idxs]

        matched_idx_map = dict(zip(initial_match_idxs, final_match_idxs))

        # Alchemical calculations require substantial portion of the two molecules
        # to be the same i.e. if a side chain is exchanged. If you are changing the whole
        # molecule, than this method is not applicable, so print a warning.

        if float(len(initial_match_idxs))/float(len(initial_mol_info.mol.GetAtoms())) < 0.3:
            print("WARNING: The size of the selected common region is less than 30% the size of the molecule at the disappearing end.")
#            sys.exit(0)

    else:

        print("No matching region could be identified!")
        sys.exit(1)

#    selected_submatch = submatches[selected_idx]
    selected_submatch = initial_match_idxs
    # Output all submatches and information informing choice
    output_MCS_file(initial_atom_info, final_atom_info,
                           matched_idx_map, output_dir)

    return (
        initial_mol_info,
        initial_atom_info,
        final_mol_info,
        final_atom_info,
        matched_idx_map,
        selected_submatch)


def correct_atom_type_mismatch(initial_match_idxs, matched_idx_map, initial_mol_info,
				 final_mol_info, initial_atom_info):
    """
    Check if the atom types for matched atoms are the same in both molecules. If not
    deletes the H-atoms (or heavy atoms + all attached Hs) with mismatches. 
    Thereafter, makes sure that the overall q_diff still falls within the specified tolerance.

    Args:
	initial_match_idxs: List of MCS atom indices for initial mol
        matched_idx_map (dict): Map index for initial molecule to index of the
                                corresponding atom in the final molecule
        initial_mol_info (MolInfo()): Structure, name and type information for
                                      the initial molecule
        final_mol_info (MolInfo()): Structure, name and type information for
                                    the final molecule
        initial_atom_info (dict): Naming and connectivity (AtomInfo) information by atom index

    Returns:
	List of indices in MCS after correction
    """

    indices_to_be_deleted = []
    new_match_idxs = []

    for init_idx, final_idx in matched_idx_map.iteritems():
	if init_idx not in indices_to_be_deleted:
	        init_type = initial_mol_info.atom_types[init_idx]
	        final_type = final_mol_info.atom_types[final_idx]
	
	        if init_type != final_type:
		    indices_to_be_deleted.append(init_idx)
		    print ("%s to be deleted from MCS due to type mismatch error." % initial_atom_info[init_idx].name)
	       	    for bound in initial_atom_info[init_idx].bound:
	                if (initial_atom_info[bound].name.strip('0123456789') == 'H') and (bound not in indices_to_be_deleted) and (bound in initial_match_idxs):
		      	   indices_to_be_deleted.append(bound)
		           print ("%s (attached to %s) to be deleted from MCS due to type mismatch error." % (initial_atom_info[bound].name, initial_atom_info[init_idx].name))
	

    new_match_idxs = list(set(initial_match_idxs) - set(indices_to_be_deleted))
    return new_match_idxs     

def delete_atoms_with_large_qdiff(initial_match_idxs, matched_idx_map,
                                initial_charge_map, final_charge_map, 
				initial_atom_info, atom_tolerance, Hatom_tolerance):
    """
    Deletes H-atoms or groups of atoms (heavy atoms + all H atoms attached to them) 
    from MCS with their qdiffs between the initial and final mols greater
    than the atom_tolerance.

    Args:
        initial_match_idxs (list): List of atom indices for MCS
        matched_idx_map (dict): Map index of atom in initial to that in the
                                final molecule
        initial_charge_map (bac_ties.bac_atom_data.PdbRdkitChargeMap): Charges for atom indices for initial molecule
        final_charge_map (bac_ties.bac_atom_data.PdbRdkitChargeMap): Charges for atom indices for final molecule
        initial_atom_info (dict): Naming and connectivity (AtomInfo) information by atom index
	atom_tolerance: Threshold value of qdiff for each heavy atom.
	Hatom_tolerance: Threshold value of qdiff for each hydrogen atom.

    Returns:
        list: List of atom indices for updated MCS with relevant atoms deleted.

    """

    indices_to_be_deleted = []
    new_match_idxs = []
    non_H_idxs = []

    for idx in initial_match_idxs:
       
       q_init = initial_charge_map.charge_rdkit_idx(idx)
       idx_final = matched_idx_map[idx]
       q_final = final_charge_map.charge_rdkit_idx(idx_final)
       q_diff = q_init - q_final
       if (initial_atom_info[idx].name.strip('0123456789') == 'H') and (abs(q_diff) > Hatom_tolerance):
	       indices_to_be_deleted.append(idx)
	       print ("%s to be deleted from MCS due to large q_diff." % initial_atom_info[idx].name)
       elif (initial_atom_info[idx].name.strip('0123456789') != 'H') and (abs(q_diff) > atom_tolerance):
	       non_H_idxs.append(idx)

    for idx in non_H_idxs:
       tmp = [idx]
       q_init = initial_charge_map.charge_rdkit_idx(idx)
       idx_final = matched_idx_map[idx]
       q_final = final_charge_map.charge_rdkit_idx(idx_final)
       for bound in initial_atom_info[idx].bound:
	    	if (initial_atom_info[bound].name.strip('0123456789') == 'H') and (bound in initial_match_idxs) and (bound not in indices_to_be_deleted):
		       tmp.append(bound)
	      	       q_init += initial_charge_map.charge_rdkit_idx(bound)
	      	       bound_final = matched_idx_map[bound]
	      	       q_final += final_charge_map.charge_rdkit_idx(bound_final)
		       
       q_diff = q_init - q_final
       if abs(q_diff) > atom_tolerance:
	   indices_to_be_deleted.extend(tmp)
	   print ("%s (and attached H-atoms) to be deleted due to large q_diff." % initial_atom_info[idx].name)
       
    
    new_match_idxs = list(set(initial_match_idxs) - set(indices_to_be_deleted))
    return new_match_idxs

def overall_qdiff_correction(initial_match_idxs, matched_idx_map, initial_charge_map, 
				final_charge_map, initial_atom_info, tolerance):
    """
    Checks if the difference between the overall charge of MCS in the
    initial and final molecules is greater than the tolerance, and if so,
    corrects it by iteratively deleting H-atoms (or groups of atoms) with largest qdiff.

    Args:
        initial_match_idxs (list): List of atom indices for MCS
        matched_idx_map (dict): Map index of atom in initial to that in the
                                final molecule
        initial_charge_map (bac_ties.bac_atom_data.PdbRdkitChargeMap): Charges for atom indices for initial molecule
        final_charge_map (bac_ties.bac_atom_data.PdbRdkitChargeMap): Charges for atom indices for final molecule
        initial_atom_info (dict): Naming and connectivity (AtomInfo) information by atom index
	tolerance: Threshold value of overall qdiff for MCS.

    Returns:
        list: List of atom indices after correction.

    """

    q_max_atom_diff = 0
    q_diff_total = 0
    
    for idx in initial_match_idxs:
       q_init = initial_charge_map.charge_rdkit_idx(idx)
       idx_final = matched_idx_map[idx]
       q_final = final_charge_map.charge_rdkit_idx(idx_final)
       q_diff = q_init - q_final
       q_diff_total += q_diff
       if abs(q_diff) > q_max_atom_diff:
	  q_max_atom_diff = abs(q_diff)
	  max_idx = idx
    
    if abs(q_diff_total) > tolerance:
       indices_to_be_deleted = []
       new_match_idxs = []
       if initial_atom_info[max_idx].name.strip('0123456789') == 'H':
	       indices_to_be_deleted.append(max_idx)
	       print ("%s to be deleted from MCS for q_total correction." % initial_atom_info[max_idx].name)
       else:
	       indices_to_be_deleted = calc_group_charge_diffs(initial_match_idxs, matched_idx_map,
						   initial_charge_map, final_charge_map, initial_atom_info)

       new_match_idxs = list(set(initial_match_idxs) - set(indices_to_be_deleted))
       tmp = []
       tmp = overall_qdiff_correction(new_match_idxs, matched_idx_map, initial_charge_map, 
					final_charge_map, initial_atom_info, tolerance)
       return tmp
    else:
       return initial_match_idxs     

def calc_group_charge_diffs(initial_match_idxs, matched_idx_map,
                                initial_charge_map, final_charge_map, initial_atom_info):
    """
    Get the atom-group (heavy atom + all H attached to it) with maximum
    charge differences between the initial and final molecules.
    Only called when the q_max_atom_diff is on a heavy atom.

    Args:
        initial_match_idxs (list): List of matched atom indices for initial mol
        matched_idx_map (dict): Map index of atom in initial to that in the
                                final molecule
        initial_charge_map (bac_ties.bac_atom_data.PdbRdkitChargeMap): Charges for atom indices for initial molecule
        final_charge_map (bac_ties.bac_atom_data.PdbRdkitChargeMap): Charges for atom indices for final molecule
        initial_atom_info (dict): Naming and connectivity (AtomInfo) information by atom index

    Returns:
        List: The atomic-group with the highest charge difference for the match

    """

    indices_to_be_deleted= []
    q_max_atom_diff = 0

    for idx in initial_match_idxs:
       if initial_atom_info[idx].name.strip('0123456789') != 'H':
	       tmp = [idx]
               q_init = initial_charge_map.charge_rdkit_idx(idx)
               idx_final = matched_idx_map[idx]
               q_final = final_charge_map.charge_rdkit_idx(idx_final)
               for bound in initial_atom_info[idx].bound:
                      if (initial_atom_info[bound].name.strip('0123456789') == 'H') and (bound in initial_match_idxs):
			tmp.append(bound)
                        q_init += initial_charge_map.charge_rdkit_idx(bound)
                        bound_final = matched_idx_map[bound]
                        q_final += final_charge_map.charge_rdkit_idx(bound_final)
 
               q_diff = q_init - q_final

               if abs(q_diff) > q_max_atom_diff:
                  q_max_atom_diff = abs(q_diff)
		  indices_to_be_deleted = tmp
		  tmp1 = idx

    print ("%s (and attached H-atoms) to be deleted from MCS for q_total correction." % initial_atom_info[tmp1].name)
    return indices_to_be_deleted

def find_disjoint_fragments(pdb_filename):
    """
    Identifies all the disjoint fragments of heavy atoms in the
    MCS and returns the atom idx for the ones to be removed.
	
    Args:
        pdb_filename: Path to a pdb file containing only the MCS

    Returns:
        RDkit Mol: The chosen (largest) disjoint fragments.

    """

    rdkit_mol = Chem.MolFromPDBFile(pdb_filename, removeHs=False)
    frag_mols = Chem.GetMolFrags(rdkit_mol, asMols=True)
    frag_idxs = Chem.GetMolFrags(rdkit_mol)
   
    if len(frag_idxs) > 1:
    	print ("\n%d disjoint fragments found in the MCS.\nChoosing the one with the largest number of heavy atoms." % len(frag_idxs))
    	frag_len_max = 0
    	max_idx = 0
    	for x in range(len(frag_idxs)):
    		if Chem.RemoveHs(frag_mols[x]).GetNumAtoms() > frag_len_max:
    			frag_len_max = Chem.RemoveHs(frag_mols[x]).GetNumAtoms()
    			max_idx = x
    
    		elif Chem.RemoveHs(frag_mols[x]).GetNumAtoms() == frag_len_max:
    			print ("Two fragments clash with the same number of heavy atoms (%d)!\nChoosing the one with higher number of rings." % Chem.RemoveHs(frag_mols[x]).GetNumAtoms())
    			mol1 = frag_mols[max_idx].GetRingInfo()
    			rings1 = mol1.NumRings()
    			mol2 = frag_mols[x].GetRingInfo()
    			rings2 = mol2.NumRings()
    			if rings2 > rings1:
    				max_idx = x
    			elif rings2 == rings1:
    				print ("The two clashing fragments also have the same number of rings!\nChoosing the one with higher number of heavy atoms across rings.\n(Preferring separated rings over fused ones.)")
    				ring_atoms1 = []
    				ring_atoms2 = []
    				for atom in mol1.AtomRings():
    					ring_atoms1.extend(x for x in atom if x not in ring_atoms1)
    				for atom in mol2.AtomRings():
    					ring_atoms2.extend(x for x in atom if x not in ring_atoms2)
    				if len(ring_atoms2) > len(ring_atoms1):
    					max_idx = x
    				elif len(ring_atoms2) == len(ring_atoms1):
    					print ("The two clashing fragments also have the same number of heavy atoms across rings!\nChoosing the one with higher number of bonds across rings.\n(Preferring separated rings over fused ones.)")
    					ring_bonds1 = []
    					ring_bonds2 = []
    					for bond in mol1.BondRings():
    						ring_bonds1.extend(x for x in bond if x not in ring_bonds1)
    					for bond in mol2.BondRings():
    						ring_bonds2.extend(x for x in bond if x not in ring_bonds2)
    					if len(ring_bonds2) > len(ring_bonds1):
    						max_idx = x
	return Chem.MolToSmarts(frag_mols[max_idx])

    else:
	return None

