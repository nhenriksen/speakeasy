import logging
import sys
import argparse

from .typing import create_unique_type_list
from .typing import write_unique_frcmod_mol2s

from .amber import write_amber_mol2
from .amber import write_amber_frcmod

from .smirnoff import mol2_to_OEMol
from .smirnoff import create_openmm_system
from .smirnoff import write_smirnoff_frcmod

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
    ap.add_argument("-i", "--input", required=True, help="Input MOL2 file with Tripos atom names")
    ap.add_argument("-fi", "--input-format", required=False,
                    default="mol2",
                    choices=["mol2", "MOL2"],
                    help="File format of input file")
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

    if args.input_format.lower() != "mol2":
        raise NotImplementedError("MOL2 is the only input format.")
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

    # Method 1
    # Go from Tripos MOL2 to GAFF MOL2 to GAFF-based FRCMOD file.
    write_amber_mol2(conversion)
    write_amber_frcmod(conversion, full=True)

    # Method 2
    # Go from Tripos MOL2 to an OpenEye OEMol
    # OpenEye OEMol to SMIRNOFF-based OpenMM System
    # OpenMM System to MOL2 and FRCMOD via ParmEd
    openeye_molecules = mol2_to_OEMol(conversion)
    topology, system = create_openmm_system(conversion, openeye_molecules)
    write_smirnoff_frcmod(conversion, topology, system)

    from .utils import map, rewrite_smirnoff_with_gaff
    atom_mapping = map(conversion)
    rewrite_smirnoff_with_gaff(conversion, atom_mapping)