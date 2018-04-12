import parmed as pmd
from openforcefield.typing.engines.smirnoff import ForceField, unit
from openeye.oechem import (
    oemolistream, oemolostream, OEIFlavor_MOL2_Forcefield,
    OEIFlavor_Generic_Default, OEIFlavor_PDB_Default, OEIFlavor_PDB_ALL,
    OEFormat_MOL2, OEFormat_MOL2H, OEWriteMolecule, OETriposAtomNames, OEMol,
    OEFormat_PDB, OESmilesToMol, OEAddExplicitHydrogens, OEHasAtomIdx,
    OEAtomGetResidue)
from aimtools.unique_types import *
import sys

class Conversion(object):
    def __init__(self):
        self.mol2_file = None
        self.output_prefix = 'mol'

    def convert(self):
        # Set OEMol
        ifs = oemolistream()
        ifs.SetFlavor(OEFormat_MOL2, OEIFlavor_MOL2_Forcefield)
        ifs.open(self.mol2_file)

        # Read in molecules
        molecules = []
        for i,mol in enumerate(ifs.GetOEMols()):
            if i > 0:
                raise Exception('Only single residue molecules are currently supported')
            OETriposAtomNames(mol)
            molecules.append(OEMol(mol))

        # Set topology
        mol2_topo = pmd.load_file(self.mol2_file, structure=True)

        # Parameterize
        ff = ForceField('forcefield/smirnoff99Frosst.offxml')
        system = ff.createSystem(
            mol2_topo.topology,
            molecules,
            nonbondedCutoff=1.1 * unit.nanometer,
            ewaldErrorTolerance=1e-4)

        # Load into Parmed
        pmd_system = pmd.openmm.topsystem.load_topology(
            mol2_topo.topology,
            system,
            mol2_topo.positions)

        parm = pmd.amber.AmberParm.from_structure(pmd_system)

        # Create unique atom types
        unique_types = create_unique_type_list(parm)

        # Write AMBER mol2 and frcmod
        write_unique_frcmod_mol2s(parm,unique_types,names=self.output_prefix)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        raise Exception('You must provide a MOL2 file as argument')

    cnvs = Conversion()
    cnvs.mol2_file = sys.argv[1]
    if len(sys.argv) == 3:
        cnvs.output_prefix = sys.argv[2]
    cnvs.convert()
