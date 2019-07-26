import os
import logging
import re
import subprocess as sp

logger = logging.getLogger(__name__)

def write_amber_mol2(conversion):
    """Use `antechamber` to convert Tripos atom types to GAFF atom types."""

    antechamber = f"""antechamber -i {conversion.input_mol2}""" +\
    f""" -fi mol2 -o {conversion.output_mol2} -fo mol2 """ +\
    f""" -at gaff2 """ +\
    f""" -dr n"""
    antechamber_output = f"{conversion.output_prefix}-ac.out"
    antechamber_input = f"{conversion.output_prefix}-ac.in"
    with open(antechamber_input, "w") as file:
        file.write("#!/usr/bin/env bash\n")
        file.write(antechamber)
    with open(antechamber_output, "w") as file:
        p = sp.Popen(["bash", antechamber_input], cwd=".", stdout=file, stderr=file)
        output, error = p.communicate()
    if p.returncode == 0:
        logger.info("MOL2 file written by antechamber.")
        # Cleanup after `antechamber`
        for temp in [
            "ANTECHAMBER_AC.AC",
            "ANTECHAMBER_AC.AC0",
            "ANTECHAMBER_BOND_TYPE.AC",
            "ANTECHAMBER_BOND_TYPE.AC0",
            "ATOMTYPE.INF",
            antechamber_input,
            antechamber_output
        ]:
            try:
                os.remove(temp)
            except OSError:
                pass
    elif p.returncode == 1:
        logger.error("Error returned by antechamber.")
        logger.error(f"Output: {output}")
        logger.error(f"Error: {error}")


def write_amber_frcmod(conversion, full=False):
    """Use `prmchk2` to to write an explicit `frcmod` file."""

    if full == False:
        p = sp.Popen(
            [
                "parmchk2",
                "-i",
                conversion.output_mol2,
                "-f",
                "mol2",
                "-o",
                conversion.output_prefix + ".frcmod",
            ]
        )
    else:
        p = sp.Popen(
            [
                "parmchk2",
                "-i",
                conversion.output_mol2,
                "-f",
                "mol2",
                "-o",
                conversion.output_prefix + ".frcmod",
                "-a",
                "Y"
            ]
        )
    output, error = p.communicate()
    if p.returncode == 0:
        logger.info("FRCMOD file written by prmchk2.")

    if p.returncode != 0:
        logger.debug("Error running `parmchk2`.")
        logger.error(f"Output: {output}")
        logger.error(f"Error: {error}")


def run_tleap(conversion):
    """Use `tleap` to create a GAFF `prmtop`"""
    tleap = f"""
    source leaprc.gaff2
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
        logger.info("`prmtop` and `inpcrd` files written by tleap.")

    if p.returncode != 0:
        logger.debug("Error running `tleap`.")
        logger.error(f"Output: {output}")
        logger.error(f"Error: {error}")
