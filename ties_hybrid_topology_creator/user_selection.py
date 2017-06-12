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
