import parmed.topologyobjects


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
            dihed_type = parmed.topologyobjects.DihedralType(phi_k, per, phase)
            typemap[key] = dihed_type
            struct.dihedral_types.append(dihed_type)
        # This is the line that has been changed (ak.bond_partners -> ai.bond_partners)
        improper = (
            ak in ai.bond_partners and aj in ai.bond_partners and al in ai.bond_partners
        )
        struct.dihedrals.append(
            parmed.topologyobjects.Dihedral(ai, aj, ak, al, improper=improper, type=dihed_type)
        )
    struct.dihedral_types.claim()
