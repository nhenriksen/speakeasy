import parmed as pmd

from openeye.oechem import (
    oemolistream, oemolostream, OEIFlavor_MOL2_Forcefield,
    OEIFlavor_Generic_Default, OEIFlavor_PDB_Default, OEIFlavor_PDB_ALL,
    OEFormat_MOL2, OEFormat_MOL2H, OEWriteMolecule, OETriposAtomNames, OEMol,
    OEFormat_PDB, OESmilesToMol, OEAddExplicitHydrogens, OEHasAtomIdx,
    OEAtomGetResidue)

from openforcefield.typing.engines.smirnoff import (
    ForceField, generateTopologyFromOEMol, generateGraphFromTopology, unit)

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
        OETriposAtomNames(mol)
        molecules.append(OEMol(mol))
    return molecules

def create_openmm_system(conversion, molecules):
    """
    Create an OpenMM system using the input MOL2 file and force field file.
    """

    topology = generateTopologyFromOEMol(molecules[0])
    ff = ForceField(conversion.ff)
    system = ff.createSystem(
        topology,
        molecules,
        nonbondedCutoff=1.1 * unit.nanometer,
        ewaldErrorTolerance=1e-4,
    )

    return topology, system

def write_smirnoff_frcmod(conversion, topology, system):
    smirnoff_prmtop = conversion.output_prefix + "-smirnoff.prmtop"
    smirnoff_frcmod = conversion.output_prefix + "-smirnoff.frcmod"
    structure = pmd.openmm.load_topology(topology, system)

    structure.save(smirnoff_prmtop)
    pmd.tools.writeFrcmod(pmd.load_file(smirnoff_prmtop),
                          smirnoff_frcmod).execute()
