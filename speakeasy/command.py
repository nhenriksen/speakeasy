import logging
import os
import sys
import argparse

from speakeasy.smirnoff import mol2_to_OEMol
from speakeasy.smirnoff import create_openmm_system, create_openmm_system_from_smiles
from speakeasy.smirnoff import write_smirnoff_frcmod

logger = logging.getLogger(__name__)

class Conversion(object):
    def __init__(self):
        self.input_mol2 = None
        self.ff = None
        self.output_mol2 = None
        self.output_frcmod = None
        self.__output_prefix = "conversion"

    @property
    def output_prefix(self):
        return ".".join(self.output_mol2.split(".")[0:-1])


def parse_arguments():
    ap = argparse.ArgumentParser(description="Speakeasy: a tool for using SMIRNOFF format force fields with Amber")
    ap.add_argument("-i", "--input", required=True, help="Input MOL2 file with Tripos atom names or SMILES string")
    ap.add_argument("-fi", "--input-format", required=False,
                    default="mol2",
                    choices=["mol2", "MOL2", "SMILES", "smiles"],
                    help="File format of input file (or 'SMILES')")
    ap.add_argument("-o", "--output", required=True, help="Output MOL2 file with Amber atom names")
    ap.add_argument("-fo", "--output-format", required=False,
                    default="mol2",
                    choices=["mol2", "MOL2"],
                    help="File format of output file")
    ap.add_argument("-s", "--smirnoff", required=True, help="SMIRNOFF format force field file location")
    ap.add_argument("-fr", "--frcmod", required=True, help="Amber FRCMOD output file name")
    ap.add_argument("--log-level",
                    default="WARNING",
                    choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                    help='Set the logging level for stderr logging',
                    )
    args = ap.parse_args()

    if args.input_format.lower() != "mol2" or args.input_format.lower() != "smiles":
        raise NotImplementedError("Only MOL2 and SMILES input are supported.")
    if args.output_format.lower() != "mol2":
        raise NotImplementedError("MOL2 is the only output format.")

    return args


def main():

    logging.captureWarnings(True)

    # Log to stderr
    logger = logging.getLogger()
    stream_handler = logging.StreamHandler(stream=sys.stderr)
    stream_handler.setFormatter(logging.Formatter('## {levelname}\n{message}', style='{'))
    logger.addHandler(stream_handler)

    args = parse_arguments()
    logger.setLevel(getattr(logging, args.log_level))

    conversion = Conversion()
    conversion.input_mol2 = args.input
    conversion.output_mol2 = args.output
    conversion.ff = args.smirnoff
    conversion.frcmod = args.frcmod

    if os.path.isfile(args.input):
        openeye_molecules = mol2_to_OEMol(conversion)
        topology, system = create_openmm_system(conversion, openeye_molecules)
    else:  # must be smiles?
        topology, system = create_openmm_system_from_smiles(conversion, args.input)

    write_smirnoff_frcmod(conversion, topology, system)