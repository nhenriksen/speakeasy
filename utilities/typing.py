import numpy as np
import parmed as pmd


def check_unique_atom_names(molecule):
    """
    Checks to make sure all the atoms in a structure have unique names. This is important for
    host molecules that are composed of multiple residues and thus, may have repeated atom names.
    Parameters
    ----------
    molecule : openeye.oechem.OEMol
    """
    logging.info('Checking all atoms have unique names...')
    atoms = molecule.GetMaxAtomIdx()
    atom_names = set()
    for atom in range(atoms):
        atom_names.add(molecule.GetAtom(OEHasAtomIdx(atom)).GetName())
    logging.debug(
        f'{atoms} atoms in structure, {len(atom_names)} unique atom names.')
    assert atoms == len(atom_names)


def _create_atom_type_list(first_chars, second_chars):
    """ Create all possible two character atom types """

    # NMH: Could make this into a list comprehension probably.

    types = []
    for first_char in first_chars:
        for second_char in second_chars:
            types.append(first_char + second_char)
    return types


def assign_unique_types(parm, equiv_ids, type_list, type_idx):
    """ Inspect the parms and assign unique atom types """

    unique_types = []

    for i, atom in enumerate(parm):

        # Determine element from atomic number
        element = str(atom.element)

        # If the equiv_id is not 0, then there is another atom to which
        # it is equivalent, so we should assign it the same uniquetype.
        if equiv_ids[i] != 0:
            # Subtract by 1 to deal with AMBER indexing --> pythong
            unique_type = unique_types[equiv_ids[i] - 1]
        else:
            # This is a new unique atom, so assign a new type
            if type_idx[element] < len(type_list[element]):
                # We still have types available that have priority names
                unique_type = type_list[element][type_idx[element]]
                type_idx[element] += 1
            else:
                # We ran out of priority names, taking from all possible.
                # First skip over any types already in use.
                while type_list['all'][type_idx['all']] in unique_types:
                    type_idx['all'] += 1
                unique_type = type_list['all'][type_idx['all']]
                type_idx['all'] += 1

        unique_types.append(unique_type)

    return unique_types


def create_unique_type_list(parms, equiv_files=None):
    """
    Create list(s) of unique atom types for each atom in the parm(s).

    Parameters:
    ----------
    parms : AmberParm or list of AmberParm
        All atoms in the AmberParm(s) will be uniquely typed.
    equiv_files : string or list of strings
        The filename(s) for the equivalent atoms file (EQUIVATOMS.DAT) that
        can be created by a patched version of antechamber. If a list, the
        ordering should match the order of the parms list. If equiv_files
        are not provided, all atoms will be considered unique. The equiv_file
        format is:

atomnum 5
atom C1 c3 1 0
atom H1 hc 2 0
atom H2 hc 3 2
atom H3 hc 4 2
atom H4 hc 5 2

        where the header (atomnum 5) gives the number of atoms, and each line
        consists of: 'atom' <atom_name> <atom_type> <index> <equiv_id>
        If the equiv_id is 0, the atom is the first of its unique type, while
        any other value indicates the serial number of the exemplar atom which
        has the same unique type.  Indexing starts at 1.

    Returns:
    -------
    unique_types : list or list of lists
        The list(s) of unique atom types for the parms
    """

    ### Maintain the convenience of not specifying as list
    if not isinstance(parms, list):
        parms = [parms]
    if not isinstance(equiv_files, list):
        equiv_files = [equiv_files]

    ### ASCII Characters to consider (don't use dash, "-"!)
    char_list = "0 1 2 3 4 5 6 7 8 9 a b c d e f g h i j k l m n o p q r\
                 s t u v w x y z A B C D E F G H I J K L M N O P Q R S T\
                 U V W X Y Z * & $ # % [ ] { } < > ? + = : ; ' . , ! ~ `\
                 @ ^ ( ) _ | / \\ \"".split()
    char_list.append(' ')

    # NMH: Are there any protected/hardcoded types that must be avoided?

    ### Generate atom type lists. We'll prioritize certain names for
    #   common elements, as determined by their atomic number.
    type_list = {}
    type_list['1'] = _create_atom_type_list(['h', 'H', '1'], char_list)
    type_list['6'] = _create_atom_type_list(['c', 'C', '6'], char_list)
    type_list['7'] = _create_atom_type_list(['n', 'N', '7'], char_list)
    type_list['8'] = _create_atom_type_list(['o', 'O', '8'], char_list)
    type_list['9'] = _create_atom_type_list(['f', 'F', '9'], char_list)
    type_list['15'] = _create_atom_type_list(['p', 'P', 'q'], char_list)
    type_list['16'] = _create_atom_type_list(['s', 'S', 'Z'], char_list)
    type_list['17'] = _create_atom_type_list(['l', 'L', '!'], char_list)
    type_list['35'] = _create_atom_type_list(['b', 'B', '$'], char_list)
    type_list['53'] = _create_atom_type_list(['i', 'I', '|'], char_list)
    # Do not use space as a first character in atom type
    type_list['all'] = _create_atom_type_list(char_list[:-1], char_list)

    ### Create a per-element index counter for atom types to bookkeep
    #   the types in use
    type_idx = {'1': 0, '6': 0, '7': 0, '8': 0, '9': 0, '15': 0, '16': 0,\
                '17': 0, '35': 0, '53': 0, 'all': 0}

    ### New unique types. Will correspond with atom index.
    unique_types = []

    for i, parm in enumerate(parms):

        # If no equiv_files are provided, assume all atoms are unique.
        # Note, the None was forced inside a list above.
        if equiv_files == [None]:
            equiv_ids = np.zeros([len(parm.atoms)], np.int32)
        else:
            equiv_file = equiv_files[i]
            # If there is not an equiv_file for this particular parm,
            # assume all atoms are unique
            if equiv_file is None:
                equiv_ids = np.zeros([len(parm.atoms)], np.int32)
            else:
                equiv_idxs,equiv_ids = np.loadtxt(equiv_file, dtype=np.int32,
                                   skiprows=1, usecols=(3,4), unpack=True)

        # Run the unique typer thing
        unique_types.append(assign_unique_types(parm, equiv_ids, type_list, type_idx))

    if len(parms) == 1:
        return unique_types[0]
    else:
        return unique_types


def write_unique_frcmod_mol2s(parms,
                              unique_types,
                              names=None,
                              batch_frcmod_file=None,
                              path='.'):
    """
    Create frcmod with terms for each unique atom.

    Parameters:
    ----------
    parms : AmberParm or list of AmberParm
        These are required to write
    unique_types :  list or list of lists
        The list(s) of unique atom types for the parms
    names : str or list of strings
        If not present, mol2/frcmods will be given default names, e.g. mol_1.mol2.
        Specifying names will custom name the mol2/frcmods.
    batch_frcmod_file : str
        If not None, the parameters for all molecules in this batch will be written
        to a single frcmod file with the specified name
    path : str
        Path to location where output files are written

    """

    # Should I add a check to make sure the parm is correct?  It needs to be
    # a structure (I think?) in order to write a mol2.

    # Allow for custom output paths
    path += '/'

    ### Maintain the convenience of not specifying as list
    if not isinstance(parms, list):
        parms = [parms]
    if not isinstance(unique_types[0], list):
        unique_types = [unique_types]
    if not isinstance(names, list):
        names = [names]

    ### Will need to keep track of frcmods that are written
    frcmod_list = []

    for i, parm in enumerate(parms):
        ### Replace atom types in raw parm data lists
        parm.parm_data['AMBER_ATOM_TYPE'] = unique_types[i]
        ### Run parmed incantation for "re-initializing" all the parameters
        parm.load_atom_info()
        parm.fill_LJ()
        ### Create parmset (organized list of the parameters)
        parmset = pmd.amber.AmberParameterSet.from_structure(parm)
        ### Write the individual molecule frcmod and mol2
        if names == [None]:
            indiv_frcmod = 'mol_' + str(i) + '.frcmod'
            mol2 = 'mol_' + str(i) + '.mol2'
        else:
            if names[i] is None:
                indiv_frcmod = 'mol_' + str(i) + '.frcmod'
                mol2 = 'mol_' + str(i) + '.mol2'
            else:
                indiv_frcmod = names[i] + '.frcmod'
                mol2 = names[i] + '.mol2'
        parmset.write(path+indiv_frcmod)
        frcmod_list.append(path+indiv_frcmod)
        parm.save(path+mol2, overwrite=True)

    ### Can't figure out a way to skip the 'write out individual frcmod
    #   and then read them all in to generate complete frcmod' approach.
    #   Don't see a way to combine parmsets.  Should look more ...
    if batch_frcmod_file:
        parmset = pmd.amber.AmberParameterSet(frcmod_list)
        parmset.write(path+batch_frcmod_file)
