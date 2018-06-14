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


def convert_to_gaff(conversion):
    # Convert SYBYL `mol2` to `GAFF` mol2
    run_antechamber(conversion)
    # Run `parmchk2` to create a `frcmod`
    run_prmchk2(conversion)
    # Run `tleap` to create a `prmtop` file with using default GAFF parameters
    run_tleap(conversion)


def run_antechamber(conversion, net_charge=0):
    "Use `antechamber` to convert SYBYL atom types to GAFF atom types."

    antechamber = f"""    
    antechamber -i {conversion.mol2_file} -fi mol2 -o {conversion.output_prefix + '-gaff.mol2'} -fo mol2 -at gaff -dr n -nc {net_charge}
    """
    antechamber_output = conversion.output_prefix + "-ac.out"
    antechamber_input = conversion.output_prefix + "-ac.in"
    with open(antechamber_input, "w") as file:
        file.write("#!/usr/bin/env bash\n")
        file.write(antechamber)
    with open(antechamber_output, "w") as file:
        p = sp.Popen(["bash", antechamber_input], cwd=".", stdout=file, stderr=file)
        output, error = p.communicate()
    if p.returncode == 0:
        logging.info("`mol2` file written by antechamber.")
        # Cleanup after `antechamber`
        for temp in [
            "ANTECHAMBER_AC.AC",
            "ANTECHAMBER_AC.AC0",
            "ANTECHAMBER_BOND_TYPE.AC",
            "ANTECHAMBER_BOND_TYPE.AC0",
            "ATOMTYPE.INF",
        ]:
            try:
                os.remove(temp)
            except OSError:
                pass
    elif p.returncode == 1:
        logging.error("Error returned by antechamber.")
        logging.error(f"Output: {output}")
        logging.error(f"Error: {error}")


def run_prmchk2(conversion):
    """Use `prmchk2` to to write an explicit `frcmod` file."""
    p = sp.Popen(
        [
            "parmchk2",
            "-i",
            conversion.mol2_file,
            "-f",
            "mol2",
            "-o",
            conversion.output_prefix + ".frcmod",
        ]
    )
    output, error = p.communicate()
    if p.returncode == 0:
        logging.info("`frcmod` file written by prmchk2.")

    if p.returncode != 0:
        logging.debug("Error running `parmchk2`.")
        logging.error(f"Output: {output}")
        logging.error(f"Error: {error}")


def run_tleap(conversion):
    """Use `tleap` to create a GAFF `prmtop`"""
    tleap = f"""
    source leaprc.gaff
    loadamberparams {conversion.output_prefix + '.frcmod'}
    mol = loadmol2 {conversion.output_prefix + '-gaff.mol2'}
    saveamberparm mol {conversion.output_prefix + '-gaff.prmtop'} {conversion.output_prefix + '-gaff.inpcrd'}
    quit
    """

    tleap_input = conversion.output_prefix + "-gaff.in"
    tleap_output = conversion.output_prefix + "-gaff.out"
    with open(tleap_input, "w") as file:
        file.write(tleap)
    with open(tleap_output, "w") as file:
        p = sp.Popen(
            ["tleap", "-f", tleap_input, ">", tleap_output],
            cwd=".",
            stdout=file,
            stderr=file,
        )
    output, error = p.communicate()
    if p.returncode == 0:
        logging.info("`prmtop` and `inpcrd` files written by tleap.")

    if p.returncode != 0:
        logging.debug("Error running `tleap`.")
        logging.error(f"Output: {output}")
        logging.error(f"Error: {error}")


def compare_parameters(reference_prmtop, target_prmtop):
    pass


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


def compare_lj_parameters(reference_prmtop, target_prmtop, reference_to_target_mapping):
    lennard_jones = pd.DataFrame()
    reference = pmd.load_file(reference_prmtop, structure=True)
    target = pmd.load_file(target_prmtop, structure=True)

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

        logging.debug(
            f"{reference_name:4} {reference_epsilon:4.3f} {reference_sigma:4.3f} → "
            f"{target_name:4} {target_epsilon:4.3f} {target_sigma:4.3f}"
        )

    return lennard_jones


def find_bonds(structure):
    df = pd.DataFrame()
    for atom in structure.atoms:
        for bond in atom.bonds:
            df = df.append(
                pd.DataFrame(
                    {
                        "atom1": bond.atom1.name,
                        "atom2": bond.atom2.name,
                        "atom1_type": bond.atom1.type,
                        "atom2_type": bond.atom2.type,
                        "req": bond.type.req,
                        "k": bond.type.k,
                    },
                    index=[0],
                ),
                ignore_index=True,
            )
    return df


def compare_bonds(reference_prmtop, target_prmtop, reference_to_target_mapping):
    pass
