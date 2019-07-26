import parmed as pmd

from openeye.oechem import (
    oemolistream, OEIFlavor_MOL2_Forcefield,
    OEFormat_MOL2, OETriposAtomNames, OEMol)

from openforcefield.typing.engines.smirnoff import ForceField
from openforcefield.topology import Molecule, Topology

from .parmed import _process_dihedral

def mol2_to_OEMol(conversion):
    """
    Convert the input MOL2 file to a list of OpenEye OEMols.
    """

    ifs = oemolistream()
    ifs.SetFlavor(OEFormat_MOL2, OEIFlavor_MOL2_Forcefield)
    ifs.open(conversion.input_mol2)

    # Read in molecules
    molecules = []
    for i, mol in enumerate(ifs.GetOEMols()):
        if i > 0:
            raise RuntimeError("Only single residue molecules are currently supported")
        molecules.append(OEMol(mol))
    return molecules

def create_openmm_system(conversion, molecules):
    """
    Create an OpenMM system using the input MOL2 file and force field file.
    """

    molecule = Molecule.from_openeye(molecules[0])
    topology = Topology.from_molecules([molecule])
    ff = ForceField(conversion.ff)
    system = ff.create_openmm_system(topology)



    return topology, system

def write_smirnoff_frcmod(conversion, topology, system, with_fixes=True):

    smirnoff_prmtop = conversion.output_prefix + "-smirnoff.prmtop"
    smirnoff_frcmod = conversion.output_prefix + "-smirnoff.frcmod"
    smirnoff_mol2 = conversion.output_prefix + "-smirnoff.mol2"


    structure = pmd.openmm.load_topology(topology.to_openmm(), system)

    pmd.openmm.topsystem._process_dihedral_old = (
        pmd.openmm.topsystem._process_dihedral
    )
    pmd.openmm.topsystem._process_dihedral = _process_dihedral


    # FIXME: ParmEd does not assign scee and scnb to the dihedrals
    # If this block is commented out, energies won't match between
    # ParmEd and tleap -generated PRMTOPS! I guess this is a bug in
    # ParmEd because the generated FRCMOD files are incorrect.
    if with_fixes:
        for dihedral in structure.dihedral_types:
            dihedral.scee = 1.2
            dihedral.scnb = 2.0

    # FIXME: OpenForceField does not provide residue names!
    structure.residues[0].name = "RES"

    # FIXME: OpenForceField does not provide atom names!
    for index, atom in enumerate(structure.atoms):
        atom.name = f"{atom.element_name}{index}"

    structure.save(smirnoff_prmtop, format="amber")
    structure.save(smirnoff_mol2)
    parameter_set = pmd.amber.AmberParameterSet.from_structure(structure)
    parameter_set.write(smirnoff_frcmod)