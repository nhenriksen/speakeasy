#!/usr/bin/env python

import pytest
from io import StringIO
from textwrap import dedent
from pathlib import Path
from subprocess import call, check_output
import shutil
import os

import numpy as np
from scipy.optimize import minimize
from openforcefield.typing.engines.smirnoff import ForceField
from openforcefield.topology import Molecule, Topology
import parmed
from parmed.topologyobjects import (
    Angle,
    AngleType,
    Atom,
    AtomType,
    Bond,
    BondType,
    Cmap,
    CmapType,
    Dihedral,
    DihedralType,
    ExtraPoint,
    Improper,
    ImproperType,
    NonbondedException,
    NonbondedExceptionType,
    RBTorsionType,
    UreyBradley,
)
from mdtraj.utils import enter_temp_directory
import sander
from simtk import openmm, unit
from simtk.openmm import app as openmm_app


TEST_FORCEFIELD = "test_forcefields/smirnoff99Frosst.offxml"


def _process_dihedral(struct, force):
    """ Adds periodic torsions to the structure """
    typemap = dict()
    for ii in range(force.getNumTorsions()):
        i, j, k, l, per, phase, phi_k = force.getTorsionParameters(ii)
        ai, aj = struct.atoms[i], struct.atoms[j]
        ak, al = struct.atoms[k], struct.atoms[l]
        key = (per, phase._value, phi_k._value)
        if key in typemap:
            dihed_type = typemap[key]
        else:
            dihed_type = DihedralType(phi_k, per, phase)
            typemap[key] = dihed_type
            struct.dihedral_types.append(dihed_type)
        # This is the line that has been changed (ak.bond_partners -> ai.bond_partners)
        improper = (
            ak in ai.bond_partners and aj in ai.bond_partners and al in ai.bond_partners
        )
        struct.dihedrals.append(
            Dihedral(ai, aj, ak, al, improper=improper, type=dihed_type)
        )
    struct.dihedral_types.claim()


############################
#
# INPUT / OUTPUT
#
############################


def smiles_to_off_objects(smiles):
    mol = Molecule.from_smiles(smiles)
    mol.generate_conformers()
    topology = Topology.from_molecules([mol])
    ff = ForceField(TEST_FORCEFIELD)
    system = ff.create_openmm_system(topology)
    return mol, topology, system


def parmed_to_prmtop(structure):
    memfile = StringIO()
    structure.save(memfile, format="amber")
    return memfile.getvalue()


def parmed_to_inpcrd(structure):
    memfile = StringIO()
    structure.save(memfile, format="rst7")
    return memfile.getvalue()


def parmed_to_frcmod(structure):
    parm_set = parmed.amber.AmberParameterSet.from_structure(structure)
    memfile = StringIO()
    parm_set.write(memfile)
    return memfile.getvalue()


def parmed_to_mol2(structure):
    memfile = StringIO()
    structure.save(memfile, format="mol2")
    return memfile.getvalue()


def tleap_to_prmtop_inpcrd(
    mol2, frcmod, forcefield="oldff/leaprc.ff99SBildn", name="mol"
):
    # source {forcefield}
    template = dedent(
        f"""
        source leaprc.gaff2
        loadamberparams {name}.frcmod
        mol = loadmol2 {name}.mol2
        saveAmberParm mol {name}.prmtop {name}.inpcrd
        quit
        """
    )
    with enter_temp_directory():
        with open("leaprc", "w") as f:
            f.write(template)
        with open(f"{name}.mol2", "w") as f:
            f.write(mol2)
        with open(f"{name}.frcmod", "w") as f:
            f.write(frcmod)
        output = check_output(["tleap"])
        with open(f"{name}.prmtop") as f:
            prmtop = f.read()
        with open(f"{name}.inpcrd") as f:
            inpcrd = f.read()
    return prmtop, inpcrd


def amber_parm(prmtop, inpcrd):
    with enter_temp_directory():
        with open("openmm.prmtop", "w") as f:
            f.write(prmtop)
        with open("openmm.inpcrd", "w") as f:
            f.write(inpcrd)
        return sander.AmberParm("openmm.prmtop", "openmm.inpcrd")


def openmm_positions(inpcrd):
    with enter_temp_directory():
        with open("openmm.inpcrd", "w") as f:
            f.write(inpcrd)
        return openmm_app.AmberInpcrdFile("openmm.inpcrd").getPositions(asNumpy=True)


def openmm_simulation(prmtop):
    with enter_temp_directory():
        with open("openmm.prmtop", "w") as f:
            f.write(prmtop)
        prmtop = openmm_app.AmberPrmtopFile("openmm.prmtop")
    system = prmtop.createSystem(
        nonbondedMethod=openmm_app.NoCutoff,
        nonbondedCutoff=99.9 * unit.nanometers,
        constraints=None,
        rigidWater=False,
    )
    integrator = openmm.VerletIntegrator(0.001)
    platform = openmm.Platform.getPlatformByName("CPU")
    return openmm_app.Simulation(prmtop.topology, system, integrator, platform)


############################
#
# Energy calculation
#
############################


def scipy_helper_sander(positions):
    sander.set_positions(positions)
    ene, frc = sander.energy_forces()
    return ene.tot * 4.184


def scipy_helper_openmm(positions, simulation):
    simulation.context.setPositions(positions.reshape(-1, 3))
    return simulation.context.getState(getEnergy=True).getPotentialEnergy()._value


def scipy_minimizer(energy_function, x0, args=()):
    return minimize(
        energy_function,
        x0,
        args,
        method="L-BFGS-B",
        jac=False,
        options=dict(maxiter=1000, disp=False, gtol=0.01),
    )


def get_energy_with_openmm(top, crd):
    simulation = openmm_simulation(top)
    positions = openmm_positions(crd)
    final = scipy_minimizer(scipy_helper_openmm, positions, simulation)
    return final["fun"]


def get_energy_with_sander(top, crd):
    inp = sander.gas_input()
    parm = amber_parm(top, crd)
    with sander.setup(parm, parm.coordinates, None, inp):
        final = scipy_minimizer(scipy_helper_sander, parm.coordinates)
        return final["fun"]


def _out_of_plane_movement(structure, distance=0.5):
    """
    Move one of the in-plane atoms out of plane. Modifies in-place!
    """
    central_atom = next(atom for atom in structure.atoms if len(atom.bond_partners) > 1)
    partner_a, partner_b = central_atom.bond_partners[:2]
    coordinates = structure.coordinates.copy()
    central_xyz = np.copy(coordinates[central_atom.idx])
    partner_a_xyz = coordinates[partner_a.idx]
    partner_b_xyz = coordinates[partner_b.idx]
    orthogonal = np.cross(partner_a_xyz - central_xyz, partner_b_xyz - central_xyz)
    normalized = orthogonal / np.linalg.norm(orthogonal)
    coordinates[central_atom.idx] += distance * normalized
    structure.coordinates = coordinates
    assert (
        np.linalg.norm(central_xyz - structure.coordinates[central_atom.idx])
        == distance
    )


@pytest.mark.parametrize("with_fixes", (True, False))
def test_benzene_improper(with_fixes):
    # In ParmEd, the central atom in an improper dihedral is the third one.
    # OpenForceField thinks it's the first one. Patch that!
    parmed.openmm.topsystem._process_dihedral_old = (
        parmed.openmm.topsystem._process_dihedral
    )
    parmed.openmm.topsystem._process_dihedral = _process_dihedral

    molecule, topology, system = smiles_to_off_objects("c1ccccc1")
    structure = parmed.openmm.load_topology(topology.to_openmm(), system)
    structure.coordinates = molecule.conformers[0]

    # Move a carbon out of the plane ring
    _out_of_plane_movement(structure, 1.0)

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

    prmtop = parmed_to_prmtop(structure)
    inpcrd = parmed_to_inpcrd(structure)

    # Compare energies between openmm and sander
    final_energy_openmm = get_energy_with_openmm(prmtop, inpcrd)
    final_energy_sander = get_energy_with_sander(prmtop, inpcrd)
    assert np.isclose(final_energy_openmm, final_energy_sander, rtol=0.001)

    # Check that tleap generates valid prmtop from frcmod
    frcmod = parmed_to_frcmod(structure)
    mol2 = parmed_to_mol2(structure)
    tleap_prmtop, tleap_inpcrd = tleap_to_prmtop_inpcrd(mol2, frcmod)

    tleap_final_energy_openmm = get_energy_with_openmm(tleap_prmtop, tleap_inpcrd)
    tleap_final_energy_sander = get_energy_with_sander(tleap_prmtop, tleap_inpcrd)
    assert np.isclose(tleap_final_energy_openmm, tleap_final_energy_sander, rtol=0.001)

    # Compare parmed and tleap prmtops with openmm
    assert np.isclose(final_energy_openmm, tleap_final_energy_openmm, rtol=0.001)
    # Compare parmed and tleap prmtops with sander
    assert np.isclose(final_energy_sander, tleap_final_energy_sander, rtol=0.001)
