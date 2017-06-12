import itertools

from rdkit.Chem import rdFMCS

from prepare_structure import *
from bac_ties.bac_atom_data import *
from user_selection import *
from output import *


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

        # This technique only works is a substantial portion of the two molecules
        # is the same i.e. if a side chain is exchanged. If you are changing the whole
        # molecule, than this method is not applicable, so the system exits.

        # At this point only a lower bound of number of matching atoms is set
        # but a more meaningful ratio could be a good alternative.

        if len(initial_match_idxs) <= 2:
            print("No sufficiently large matching region (excluding incomplete rings) could be identified!")
            sys.exit(0)

        check_atom_type_match(
            matched_idx_map,
            initial_mol_info,
            final_mol_info)

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
            filtered_q_total = [
                x for x in all_options if q_diffs[x] < tolerance]
            filtered_q_atom = [
                x for x in filtered_q_total if q_max_atom_diffs[x] < atom_tolerance]

            all_filtered_matches = []

            if len(filtered_q_atom) > 0:

                for submatch_idx in filtered_q_atom:

                    submatch = submatches[submatch_idx]

                    no_heavy = 0

                    for idx in submatch:

                        if initial_atom_info[idx].name.strip(
                                '0123456789') != 'H':
                            no_heavy += 1

                    if no_heavy > 1:
                        all_filtered_matches.append(submatch_idx)

                if len(all_filtered_matches) > 0:

                    selected_idx = all_filtered_matches[0]

                else:

                    print(
                        'Warning: Only acceptable matches with at most 1 heavy atom found')
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

    return (
        initial_mol_info,
        initial_atom_info,
        final_mol_info,
        final_atom_info,
        matched_idx_map,
        selected_submatch)


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

                bound_match = list(
                    set(info.bound).intersection(set(updated_match)))

                for bound_idx in bound_match:

                    for check_idx in range(len(disjoint_sections)):

                        if bound_idx in disjoint_sections[check_idx]:
                            section_idx = check_idx
                            break

                # Add this atom and any bound to it to appropriate list
                if section_idx is None:

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

            print(
                'WARNING: Matched atoms of different types: {0:s} ({1:s}) - {2:s} ({3:s})'.format(
                    init_name,
                    init_type,
                    final_name,
                    final_type))
            type_match_error = True

    if type_match_error:

        choice = input(
            "Given type match error do you want to proceed (y/[n])?")
        if choice.lower() != 'y':
            sys.exit(1)

    return


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

                    potential_stop_idxs = [
                        x for x in bridge_info.bound if x in match_idxs]

                # Get selections of atoms to remove for each potential section
                # to be removed
                for stop_idx in potential_stop_idxs:

                    stop_idxs = [stop_idx]

                    remove_selection = sorted(get_linked_to_remove(bridge_idx,
                                                                   atom_info,
                                                                   match_idxs,
                                                                   stop_idxs))

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
                bridge_idxs += [
                    x for x in get_bridge_atoms(
                        mol, new_match) if x not in trial_bridges]

            if len(set(bridge_idxs)) < 1:
                break

    # Need a comprehensive list of all possible removals that can be used to
    # produce a valid submatch
    # a) Create all possible combinations of section removals
    tmp = []

    for comb_length in range(1, len(options) + 1):

        for combination in itertools.combinations(options, comb_length):

            for product in itertools.product(*combination):

                combined_removal_atoms = set(
                    [item for sublist in product for item in sublist])

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

        submatch = sorted(set(match_idxs) - set(removal_list))

        if submatch:
            submatches.append(submatch)

    # for i in range(len(submatches)):
    #     print(submatches[i])
    #     print(full_options[i])

    # sys.exit()

    return submatches


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


def get_bridge_ring(ring_no, atom_rings, atom_info):
    idxs_off_ring = []

    seen_rings = [ring_no]
    check_ring_atoms = atom_rings[ring_no]

    while len(check_ring_atoms):

        check_rings = []

        for idx in check_ring_atoms:

            bound = atom_info[idx].bound

            for bound_idx in bound:

                if bound_idx not in check_ring_atoms:

                    if atom_info[bound_idx].rings:

                        for new_ring in atom_info[idx].rings:

                            if new_ring not in check_rings and new_ring not in seen_rings:
                                check_rings.append(new_ring)

                    else:

                        idxs_off_ring.append(bound_idx)

        new_check_atoms = []

        for ring_idx in check_rings:
            for idx in atom_rings[ring_idx]:
                if idx not in check_ring_atoms:
                    new_check_atoms.append(idx)

        check_ring_atoms = new_check_atoms

        seen_rings += check_rings

    return idxs_off_ring


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

    start_ring = atom_info[idx].rings[0]

    off_idxs = get_bridge_ring(start_ring, atom_rings, atom_info)

    for off_idx in off_idxs:
        if off_idx in matched:
            stop_idxs.append(off_idx)

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

