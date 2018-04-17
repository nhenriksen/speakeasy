import parmed as pmd
from openforcefield.typing.engines.smirnoff import ForceField, unit
from openeye.oechem import (
    oemolistream, oemolostream, OEIFlavor_MOL2_Forcefield,
    OEIFlavor_Generic_Default, OEIFlavor_PDB_Default, OEIFlavor_PDB_ALL,
    OEFormat_MOL2, OEFormat_MOL2H, OEWriteMolecule, OETriposAtomNames, OEMol,
    OEFormat_PDB, OESmilesToMol, OEAddExplicitHydrogens, OEHasAtomIdx,
    OEAtomGetResidue)
import aimtools 
import sys

class Conversion(object):
    def __init__(self):
        self.mol2_file = None
        self.output_prefix = 'mol'
        self.mol2_topo = None
        self.off_system = None
        self.pmd_system = None
        self.parm = None
        self.molecules = []
        self.labels = None

        

    def convert(self):
        # Set OEMol
        ifs = oemolistream()
        ifs.SetFlavor(OEFormat_MOL2, OEIFlavor_MOL2_Forcefield)
        ifs.open(self.mol2_file)

        # Read in molecules
        for i,mol in enumerate(ifs.GetOEMols()):
            if i > 0:
                raise Exception('Only single residue molecules are currently supported')
            OETriposAtomNames(mol)
            self.molecules.append(OEMol(mol))

        # Set topology
        self.mol2_topo = pmd.load_file(self.mol2_file, structure=True)

        # Parameterize
        ff = ForceField('forcefield/smirnoff99Frosst.offxml')
        self.labels = ff.labelMolecules(self.molecules, verbose=False)
        self.off_system = ff.createSystem(
            self.mol2_topo.topology,
            self.molecules,
            nonbondedCutoff=1.1 * unit.nanometer,
            ewaldErrorTolerance=1e-4)

        # Load into Parmed
        self.pmd_system = pmd.openmm.topsystem.load_topology(
            self.mol2_topo.topology,
            self.off_system,
            self.mol2_topo.positions)

        # Convert to AmberParm
        self.parm = pmd.amber.AmberParm.from_structure(self.pmd_system)

        # HACKY PART!!
        # Amber specifies that the third atom in an improper is the central
        # atom, but smirnoff currently specifies the second atom. A check for
        # impropers was conducted during pmd.openmm.topsystem.load_topology(),
        # but that looked at the third atom, so we'll recheck the second atom.
        for i,dihedral in enumerate(cnvs.parm.dihedrals):
            a1 = dihedral.atom1
            a2 = dihedral.atom2
            a3 = dihedral.atom3
            a4 = dihedral.atom4
            if a1 in a2.bond_partners and a3 in a2.bond_partners and a4 in a2.bond_partners:
                (dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4) = (a3,a4,a2,a1)
                dihedral.improper = True

        # Create unique atom types
        unique_types = aimtools.unique_types.create_unique_type_list(self.parm)

        # Write AMBER mol2 and frcmod
        aimtools.unique_types.write_unique_frcmod_mol2s(self.parm,unique_types,names=self.output_prefix)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        raise Exception('You must provide a MOL2 file as argument')

    cnvs = Conversion()
    cnvs.mol2_file = sys.argv[1]
    if len(sys.argv) == 3:
        cnvs.output_prefix = sys.argv[2]
    cnvs.convert()
