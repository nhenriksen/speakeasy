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

    antechamber = f"""antechamber -i {conversion.mol2_file}"""
    f""" -fi mol2 -o {conversion.output_prefix + '-gaff.mol2'} -fo mol2 """
    f""" -at gaff """
    f""" -dr n -nc {net_charge}"""
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
