import logging
import re
from .smirnoff import mol2_to_OEMol

from openeye.oechem import (
    oemolistream, oemolostream, OEIFlavor_MOL2_Forcefield,
    OEIFlavor_Generic_Default, OEIFlavor_PDB_Default, OEIFlavor_PDB_ALL,
    OEFormat_MOL2, OEFormat_MOL2H, OEWriteMolecule, OETriposAtomNames, OEMol,
    OEFormat_PDB, OESmilesToMol, OEAddExplicitHydrogens, OEHasAtomIdx,
    OEAtomGetResidue)

from openforcefield.typing.engines.smirnoff import (
    generateTopologyFromOEMol, generateGraphFromTopology)
from networkx.algorithms import isomorphism

logger = logging.getLogger(__name__)

def mol2_to_OEMol(file):
    """
    Convert the input MOL2 file to a list of OpenEye OEMols.
    """

    ifs = oemolistream()
    ifs.SetFlavor(OEFormat_MOL2, OEIFlavor_MOL2_Forcefield)
    ifs.open(file)

    # Read in molecules
    molecules = []
    for i, mol in enumerate(ifs.GetOEMols()):
        if i > 0:
            raise RuntimeError("Only single residue molecules are currently supported")
        OETriposAtomNames(mol)
        molecules.append(OEMol(mol))
    return molecules


def map(conversion):

    logger.debug(f'Generating map between atoms...')

    # GAFF MOL2
    reference_mol = mol2_to_OEMol(conversion.output_mol2)[0]
    # SMIRNOFF MOL2
    target_mol = mol2_to_OEMol(conversion.output_prefix + "-smirnoff.mol2")[0]

    reference_topology = generateTopologyFromOEMol(reference_mol)
    target_topology = generateTopologyFromOEMol(target_mol)

    reference_graph = generateGraphFromTopology(reference_topology)
    target_graph = generateGraphFromTopology(target_topology)

    reference_to_target_mapping = dict()
    target_to_reference_mapping = dict()
    graph_matcher = isomorphism.GraphMatcher(reference_graph, target_graph)
    if graph_matcher.is_isomorphic():
        logging.debug('Determining mapping...')
        logging.debug('Reference → Target')
        for (reference_atom, target_atom) in graph_matcher.mapping.items():
            reference_to_target_mapping[reference_atom] = target_atom
            reference_name = reference_mol.GetAtom(
                OEHasAtomIdx(reference_atom)).GetName()
            target_name = target_mol.GetAtom(
                OEHasAtomIdx(target_atom)).GetName()
            reference_type = reference_mol.GetAtom(
                OEHasAtomIdx(reference_atom)).GetType()
            target_type = target_mol.GetAtom(
                OEHasAtomIdx(target_atom)).GetType()

            target_to_reference_mapping[target_type] = reference_type
            logging.debug(f'{reference_name:5} {reference_type:3} {reference_atom:3d} → '
                          f'{target_atom:3d} {target_type:3} {target_name:5}')
    else:
        logging.error('Graph is not isomorphic.')

    return target_to_reference_mapping


def remap_atoms(smirnoff_frcmod, target_to_reference_mapping):
    opening = "MASS"
    between = "(.*?)"
    closing = "BOND"
    pattern = opening + between + closing
    atoms = re.findall(pattern, smirnoff_frcmod, re.DOTALL)
    atom_list = atoms[0].split("\n")
    new_atoms = []
    for atom in atom_list[1:-2]:
        atom_type, mass = atom.split()
        new_atoms.append(f"{target_to_reference_mapping[atom_type]}{mass:>10}")

    return new_atoms

def remap_bonds(smirnoff_frcmod, target_to_reference_mapping):
    opening = "BOND"
    between = "(.*?)"
    closing = "ANGLE"
    pattern = opening + between + closing
    bonds = re.findall(pattern, smirnoff_frcmod, re.DOTALL)
    bond_list = bonds[0].split("\n")
    new_bonds = []
    for bond in bond_list[1:-2]:
        atom_one, atom_two, k, r = bond.split()
        atom_two = atom_two.replace("-", "")
        new_atom_one = target_to_reference_mapping[atom_one]
        new_atom_two = target_to_reference_mapping[atom_two]
        new_bonds.append(f"{new_atom_one:<2}-{new_atom_two:<2}{k:>12}{r:>7}")

    return new_bonds

def remap_angles(smirnoff_frcmod, target_to_reference_mapping):
    opening = "ANGLE"
    between = "(.*?)"
    closing = "DIHE"
    pattern = opening + between + closing
    angles = re.findall(pattern, smirnoff_frcmod, re.DOTALL)
    angles_list = angles[0].split("\n")

    new_angles = []
    for angle in angles_list[1:-2]:
        atom_one, atom_two, atom_three, k, theta = angle.split()
        atom_two = atom_two.replace("-", "")
        atom_three = atom_three.replace("-", "")
        new_atom_one = target_to_reference_mapping[atom_one]
        new_atom_two = target_to_reference_mapping[atom_two]
        new_atom_three = target_to_reference_mapping[atom_three]
        new_angles.append(f"{new_atom_one:<2}-{new_atom_two:<2}-{new_atom_three:<2}{k:>11}{theta:>9}")

    return new_angles

def remap_dihedrals(smirnoff_frcmod, target_to_reference_mapping):
    opening = "DIHE"
    between = "(.*?)"
    closing = "IMPROPER"
    pattern = opening + between + closing
    dihedrals = re.findall(pattern, smirnoff_frcmod, re.DOTALL)
    dihedrals_list = dihedrals[0].split("\n")
    new_dihedrals = []
    for dihedral in dihedrals_list[1:-2]:
        try:
            atom_one, atom_two, atom_three, atom_four, idiv, k, theta, per, scee, scnb = dihedral.split()
        except ValueError as e:
            if str(e) == "not enough values to unpack (expected 10, got 9)":
                atom_one, atom_two, atom_three, stuff = dihedral.split("-")
                atom_four, idiv, k, theta, per, scee, scnb = stuff.split()


        atom_two = atom_two.replace("-", "")
        atom_three = atom_three.replace("-", "")
        atom_four = atom_four.replace("-", "")
        new_atom_one = target_to_reference_mapping[atom_one.strip()]
        new_atom_two = target_to_reference_mapping[atom_two.strip()]
        new_atom_three = target_to_reference_mapping[atom_three.strip()]
        new_atom_four = target_to_reference_mapping[atom_four.strip()]
        new_dihedrals.append(
            f"{new_atom_one:<2}-{new_atom_two:<2}-{new_atom_three:<2}-{new_atom_four:<2}{idiv:>5}{k:>15}{theta:>9}{per:>6}    {scee:<6} {scnb}")

    return new_dihedrals

def remap_impropers(smirnoff_frcmod, target_to_reference_mapping):
    raise NotImplementedError

def remap_nonb(smirnoff_frcmod,target_to_reference_mapping):
    opening = "NONB"
    forwards = ".*$"
    pattern = opening + forwards
    nonb = re.findall(pattern, smirnoff_frcmod, re.DOTALL)
    nonb_list = nonb[0].split("\n")
    new_nonb = []
    for nonb in nonb_list[1:-2]:
        atom_type, sigma, epsilon = nonb.split()
        new_nonb.append(f"{target_to_reference_mapping[atom_type]}{sigma:>14}{epsilon:>14}")

    return new_nonb

def rewrite_smirnoff_with_gaff(conversion, target_to_reference_mapping):
    with open(conversion.output_prefix + "-smirnoff.frcmod", "r") as file:
        smirnoff_frcmod = file.read()

    atoms = remap_atoms(smirnoff_frcmod, target_to_reference_mapping)
    bonds = remap_bonds(smirnoff_frcmod, target_to_reference_mapping)
    angles = remap_angles(smirnoff_frcmod, target_to_reference_mapping)
    dihedrals = remap_dihedrals(smirnoff_frcmod,target_to_reference_mapping)
    # Impropers not implemented yet...
    nonb = remap_nonb(smirnoff_frcmod, target_to_reference_mapping)

    with open(conversion.output_prefix + "-smirnoff-remapped.frcmod", "w") as file:
        file.write(smirnoff_frcmod.split("\n")[0])
        file.write("\n")
        file.write("MASS\n")
        for atom in atoms:
            file.write(atom)
            file.write("\n")
        file.write("\n")
        file.write("BOND\n")
        for bond in bonds:
            file.write(bond)
            file.write("\n")
        file.write("\n")
        file.write("ANGLE\n")
        for angle in angles:
            file.write(angle)
            file.write("\n")
        file.write("\n")
        file.write("DIHE\n")
        for dihedral in dihedrals:
            file.write(dihedral)
            file.write("\n")
        file.write("\n")
        file.write("IMPROPER\n")
        file.write("\n")
        # Impropers
        file.write("NONB\n")
        for atom in nonb:
            file.write(atom)
            file.write("\n")
        file.write("\n")

