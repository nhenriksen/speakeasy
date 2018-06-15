import numpy as np
import os as os
import logging as logging
import pandas as pd
import subprocess as sp
import parmed as pmd
from networkx.algorithms import isomorphism
from openforcefield.typing.engines.smirnoff import (
    ForceField,
    generateTopologyFromOEMol,
    generateGraphFromTopology,
)
from openeye.oechem import (
    oemolistream,
    oemolostream,
    OEIFlavor_MOL2_Forcefield,
    OEIFlavor_Generic_Default,
    OEIFlavor_PDB_Default,
    OEIFlavor_PDB_ALL,
    OEFormat_MOL2,
    OEFormat_MOL2H,
    OEWriteMolecule,
    OETriposAtomNames,
    OEMol,
    OEFormat_PDB,
    OESmilesToMol,
    OEAddExplicitHydrogens,
    OEHasAtomIdx,
    OEAtomGetResidue,
)


def compare_parameters(reference_prmtop, target_prmtop):
    print("Establishing mapping between structures...")
    reference_to_target_mapping = create_atom_map(reference_prmtop, target_prmtop)
    reference = pmd.load_file(reference_prmtop, structure=True)
    target = pmd.load_file(target_prmtop, structure=True)

    print("Comparing LJ parameters...")
    lj = compare_lj_parameters(reference, target, reference_to_target_mapping)
    print("Comparing bond parameters...")
    bonds = compare_bonds(reference, target, reference_to_target_mapping)
    print("Comparing angle parameters...")
    angles = compare_angles(reference, target, reference_to_target_mapping)
    print("Comparing dihedral parameters...")
    dihedrals = compare_dihedrals(reference, target, reference_to_target_mapping)
    # print("Comparing improper parameters...")


def create_atom_map(reference_prmtop, target_prmtop):
    reference = pmd.load_file(reference_prmtop, structure=True)
    target = pmd.load_file(target_prmtop, structure=True)

    reference_graph = generateGraphFromTopology(reference.topology)
    target_graph = generateGraphFromTopology(target.topology)

    graph_matcher = isomorphism.GraphMatcher(reference_graph, target_graph)

    reference_to_target_mapping = dict()

    if graph_matcher.is_isomorphic():
        logging.debug("Reference → Target")
        for (reference_atom, target_atom) in graph_matcher.mapping.items():

            reference_to_target_mapping[reference_atom] = target_atom

            reference_name = reference[reference_atom].name
            target_name = target[target_atom].name

            reference_type = reference[reference_atom].type
            target_type = target[target_atom].type

            # ParmEd is 0-indexed.
            # Add 1 to match AMBER-style indexing.
            logging.debug(
                f"{reference_name:4} {reference_type:4} {reference_atom + 1:3d} → "
                f"{target_atom + 1:3d} {target_type:4} {target_name:4}"
            )

    return reference_to_target_mapping


def find_bonds(structure):
    df = pd.DataFrame()
    for atom in structure.atoms:
        for bond in atom.bonds:
            df = df.append(
                pd.DataFrame(
                    {
                        "atom1": bond.atom1.name,
                        "atom2": bond.atom2.name,
                        "atom1_idx": bond.atom1.idx,
                        "atom2_idx": bond.atom2.idx,
                        "atom1_type": bond.atom1.type,
                        "atom2_type": bond.atom2.type,
                        "req": bond.type.req,
                        "k": bond.type.k,
                    },
                    index=[0],
                ),
                ignore_index=True,
            )
    return df.drop_duplicates()


def find_angles(structure):
    df = pd.DataFrame()
    for atom in structure.atoms:
        for angle in atom.angles:
            df = df.append(
                pd.DataFrame(
                    {
                        "atom1": angle.atom1.name,
                        "atom2": angle.atom2.name,
                        "atom3": angle.atom3.name,
                        "atom1_idx": angle.atom1.idx,
                        "atom2_idx": angle.atom2.idx,
                        "atom3_idx": angle.atom3.idx,
                        "atom1_type": angle.atom1.type,
                        "atom2_type": angle.atom2.type,
                        "atom3_type": angle.atom3.type,
                        "thetaeq": angle.type.theteq,
                        "k": angle.type.k,
                    },
                    index=[0],
                ),
                ignore_index=True,
            )
    return df.drop_duplicates()


def find_dihedrals(structure):
    df = pd.DataFrame()
    for atom in structure.atoms:
        for dihedral in atom.dihedrals:
            df = df.append(
                pd.DataFrame(
                    {
                        "atom1": dihedral.atom1.name,
                        "atom2": dihedral.atom2.name,
                        "atom3": dihedral.atom3.name,
                        "atom4": dihedral.atom4.name,
                        "atom1_idx": dihedral.atom1.idx,
                        "atom2_idx": dihedral.atom2.idx,
                        "atom3_idx": dihedral.atom3.idx,
                        "atom4_idx": dihedral.atom4.idx,
                        "atom1_type": dihedral.atom1.type,
                        "atom2_type": dihedral.atom2.type,
                        "atom3_type": dihedral.atom3.type,
                        "atom4_type": dihedral.atom4.type,
                        "per": dihedral.type.per,
                        "phi_k": dihedral.type.phi_k,
                        "phase": dihedral.type.phase,
                    },
                    index=[0],
                ),
                ignore_index=True,
            )
    return df.drop_duplicates()


def label_smirks(structure_mol2):
    ifs = oemolistream()
    flavor = OEIFlavor_MOL2_Forcefield
    ifs.SetFlavor(OEFormat_MOL2, flavor)

    molecules = []
    # Read in molecules
    for mol in ifs.GetOEMols():
        OETriposAtomNames(mol)
        # Add all the molecules in this file to a list, but only return the first one.
        molecules.append(OEMol(mol))
        # This should now handle single-residue and multi-residue hosts.

    # Parameterize
    ff = ForceField("forcefield/smirnoff99Frosst.offxml")
    labels = ff.labelMolecules(molecules, verbose=True)
    return molecules, labels


def compare_lj_parameters(reference, target, reference_to_target_mapping, verbose=True):
    lennard_jones = pd.DataFrame()

    logging.debug("Reference → Target")
    logging.debug(
        f"{'Name':4} {'Eps':5} {'Sigma':5} → " f"{'Name':4} {'Eps':5} {'Sigma':5}"
    )
    for reference_atom, target_atom in reference_to_target_mapping.items():

        reference_name = reference[reference_atom].name
        reference_type = reference[reference_atom].type
        reference_sigma = reference[reference_atom].sigma
        reference_epsilon = reference[reference_atom].epsilon

        target_name = target[target_atom].name
        target_type = target[target_atom].type
        target_sigma = target[target_atom].sigma
        target_epsilon = target[target_atom].epsilon

        lennard_jones = lennard_jones.append(
            pd.DataFrame(
                {
                    "target_name": target_name,
                    "target_type": target_type,
                    "target_e": np.round(target_epsilon, decimals=5),
                    "target_s": np.round(target_sigma, decimals=5),
                    "reference_name": reference_name,
                    "reference_type": reference_type,
                    "reference_e": np.round(reference_epsilon, decimals=5),
                    "reference_s": np.round(reference_sigma, decimals=5),
                },
                index=[0],
            ),
            ignore_index=True,
        )

        if verbose:
            if (np.round(reference_epsilon, 4) != np.round(target_epsilon, 4)) or (
                np.round(reference_sigma, 4) != np.round(target_sigma, 4)
            ):
                print(
                    f"\x1b[31m{reference_name:>4} {reference_sigma:4.3f} {reference_epsilon:4.3f} → "
                    f"{target_name:>4} {target_sigma:4.3f} {target_epsilon:4.3f}\x1b[0m"
                )
            else:
                print(
                    f"{reference_name:>4} {reference_sigma:4.3f} {reference_epsilon:4.3f} → "
                    f"{target_name:>4} {target_sigma:4.3f} {target_epsilon:4.3f}"
                )

    return lennard_jones


def compare_bonds(reference, target, reference_to_target_mapping, verbose=True):
    reference_bonds = find_bonds(reference)
    target_bonds = find_bonds(target)

    assert len(reference.bonds) == len(reference_bonds)
    assert len(target.bonds) == len(target_bonds)

    bonds = pd.DataFrame()
    for reference_atom, target_atom in reference_to_target_mapping.items():

        reference_atom_bonds = reference_bonds[
            (reference_bonds["atom1_idx"] == reference_atom)
            | (reference_bonds["atom2_idx"] == reference_atom)
        ]
        target_atom_bonds = target_bonds[
            (target_bonds["atom1_idx"] == target_atom)
            | (target_bonds["atom2_idx"] == target_atom)
        ]
        df = reference_atom_bonds.join(target_atom_bonds, lsuffix="_r", rsuffix="_t")
        for index, bond in df.iterrows():

            reference_atom1 = bond["atom1_r"]
            reference_atom2 = bond["atom2_r"]
            reference_k = bond["k_r"]
            reference_req = bond["req_r"]

            target_atom1 = bond["atom1_t"]
            target_atom2 = bond["atom2_t"]
            target_k = bond["k_t"]
            target_req = bond["req_t"]

            if verbose:
                if (np.round(reference_k, 4) != np.round(target_k, 4)) or (
                    np.round(reference_req, 4) != np.round(target_req, 4)
                ):
                    print(
                        f"\x1b[31m{reference_atom1:>4}--{reference_atom2:<4} {reference_k:4.3f} {reference_req:4.3f} → "
                        f"{target_atom1:>4}--{target_atom2:<4} {target_k:4.3f} {target_req:4.3f}\x1b[0m"
                    )
                else:
                    print(
                        f"{reference_atom1:>4}--{reference_atom2:<4} {reference_k:4.3f} {reference_req:4.3f} → "
                        f"{target_atom1:>4}--{target_atom2:<4} {target_k:4.3f} {target_req:4.3f}"
                    )
        bonds = bonds.append(df, ignore_index=True)

    return bonds


def compare_angles(reference, target, reference_to_target_mapping, verbose=True):
    reference_angles = find_angles(reference)
    target_angles = find_angles(target)

    assert len(reference.angles) == len(reference_angles)
    assert len(target.angles) == len(target_angles)

    angles = pd.DataFrame()
    for reference_atom, target_atom in reference_to_target_mapping.items():

        reference_atom_angles = reference_angles[
            (reference_angles["atom1_idx"] == reference_atom)
            | (reference_angles["atom2_idx"] == reference_atom)
            | (reference_angles["atom3_idx"] == reference_atom)
        ]
        target_atom_angles = target_angles[
            (target_angles["atom1_idx"] == target_atom)
            | (target_angles["atom2_idx"] == target_atom)
            | (target_angles["atom3_idx"] == target_atom)
        ]
        df = reference_atom_angles.join(target_atom_angles, lsuffix="_r", rsuffix="_t")
        for index, angle in df.iterrows():

            reference_atom1 = angle["atom1_r"]
            reference_atom2 = angle["atom2_r"]
            reference_atom3 = angle["atom3_r"]
            reference_k = angle["k_r"]
            reference_thetaeq = angle["thetaeq_r"]

            target_atom1 = angle["atom1_t"]
            target_atom2 = angle["atom2_t"]
            target_atom3 = angle["atom3_t"]
            target_k = angle["k_t"]
            target_thetaeq = angle["thetaeq_t"]

            if verbose:
                if (np.round(reference_k, 4) != np.round(target_k, 4)) or (
                    np.round(target_thetaeq, 4) != np.round(target_thetaeq, 4)
                ):
                    print(
                        f"\x1b[31m{reference_atom1:>4}--{reference_atom2:<4}--{reference_atom3:<4} {reference_k:4.3f} {reference_thetaeq:4.3f} → "
                        f"{target_atom1:>4}--{target_atom2:<4}--{target_atom3:<4} {target_k:4.3f} {target_thetaeq:4.3f}\x1b[0m"
                    )
                else:
                    print(
                        f"{reference_atom1:>4}--{reference_atom2:<4}--{reference_atom3:<4} {reference_k:4.3f} {reference_thetaeq:4.3f} → "
                        f"{target_atom1:>4}--{target_atom2:<4}--{target_atom3:<4} {target_k:4.3f} {target_thetaeq:4.3f}"
                    )
        angles = angles.append(df, ignore_index=True)

    return angles


def compare_dihedrals(reference, target, reference_to_target_mapping, verbose=True):
    reference_dihedrals = find_dihedrals(reference)
    target_dihedrals = find_dihedrals(target)

    assert len(reference.dihedrals) == len(reference_dihedrals)
    assert len(target.dihedrals) == len(target_dihedrals)

    dihedrals = pd.DataFrame()
    for reference_atom, target_atom in reference_to_target_mapping.items():

        reference_atom_dihedrals = reference_dihedrals[
            (reference_dihedrals["atom1_idx"] == reference_atom)
            | (reference_dihedrals["atom2_idx"] == reference_atom)
            | (reference_dihedrals["atom3_idx"] == reference_atom)
            | (reference_dihedrals["atom4_idx"] == reference_atom)
        ]
        target_atom_dihedrals = target_dihedrals[
            (target_dihedrals["atom1_idx"] == target_atom)
            | (target_dihedrals["atom2_idx"] == target_atom)
            | (target_dihedrals["atom3_idx"] == target_atom)
            | (target_dihedrals["atom4_idx"] == target_atom)
        ]

        # Check which atom matches and then join the DataFrames on that atom and periodicity.
        reference_match = None
        for atom in ["atom1_idx", "atom2_idx", "atom3_idx", "atom4_idx"]:
            if (reference_atom_dihedrals[atom] == reference_atom).any():
                reference_match = atom
            else:
                pass

        target_match = None
        for atom in ["atom1_idx", "atom2_idx", "atom3_idx", "atom4_idx"]:
            if (target_atom_dihedrals[atom] == target_atom).any():
                target_match = atom
            else:
                pass
        # Check that the atom matches in the same location
        assert reference_match == target_match

        df = reference_atom_dihedrals.merge(
            target_atom_dihedrals,
            suffixes=("_r", "_t"),
            left_on="per",
            right_on="per",
            # left_on=[reference_match, "per"],
            # right_on=[target_match, "per"],
        )
        for index, dihedral in df.iterrows():

            # Ensure we are comparing the same atoms.
            # We are not handling the case where the order of the atoms in the dihedral is reversed.

            atoms = ["atom1_idx", "atom2_idx", "atom3_idx", "atom4_idx"]
            atoms_r = [
                reference_to_target_mapping[dihedral[f"{atom}_r"]] for atom in atoms
            ]
            atoms_t = [dihedral[f"{atom}_t"] for atom in atoms]
            if not all([i == j for i, j in zip(atoms_r, atoms_t)]):
                continue

            reference_atom1 = dihedral["atom1_r"]
            reference_atom2 = dihedral["atom2_r"]
            reference_atom3 = dihedral["atom3_r"]
            reference_atom4 = dihedral["atom4_r"]
            reference_k = dihedral["phi_k_r"]
            reference_per = dihedral["per"]
            reference_phase = dihedral["phase_r"]

            target_atom1 = dihedral["atom1_t"]
            target_atom2 = dihedral["atom2_t"]
            target_atom3 = dihedral["atom3_t"]
            target_atom4 = dihedral["atom4_t"]
            target_k = dihedral["phi_k_t"]
            target_per = dihedral["per"]
            target_phase = dihedral["phase_t"]

            if verbose:
                if (np.round(reference_k, 3) != np.round(target_k, 3)) or (
                    np.round(target_phase, 3) != np.round(target_phase, 3)
                ):
                    print(
                        f"\x1b[31m{reference_atom1:>4}--{reference_atom2:<4}--{reference_atom3:<4}--{reference_atom4:<4} {reference_k:4.3f} {reference_phase:4.0f} {reference_per:4.0f} → "
                        f"{target_atom1:>4}--{target_atom2:<4}--{target_atom3:<4}--{target_atom4:<4} {target_k:4.3f} {target_phase:4.0f} {target_per:4.0f}\x1b[0m"
                    )
                else:
                    print(
                        f"{reference_atom1:>4}--{reference_atom2:<4}--{reference_atom3:<4}--{reference_atom4:<4} {reference_k:4.3f} {reference_phase:4.0f} {reference_per:4.0f} → "
                        f"{target_atom1:>4}--{target_atom2:<4}--{target_atom3:<4}--{target_atom4:<4} {target_k:4.3f} {target_phase:4.0f} {target_per:4.0f} "
                    )
        dihedrals = dihedrals.append(df, ignore_index=True)

    return dihedrals
