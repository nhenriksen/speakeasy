def get_net_charge(conversion):
    """Extract the net charge from the provided `mol2` file."""

    opening = re.escape(r"@<TRIPOS>ATOM")
    between = "(.*?)"
    closing = re.escape(r"@<TRIPOS>BOND")
    pattern = opening + between + closing

    with open(conversion.input_mol2, "r") as file:
        data = file.read()
    atoms = re.findall(pattern, data, re.DOTALL)
    if len(atoms) < 1:
        logger.error(f"Could not extract atoms from {conversion.input_mol2}")
    atom_list = atoms[0].split("\n")[1:-1]
    net_charge = 0
    for atom in atom_list:
        fields = atom.split()
        net_charge += float(fields[-1])
    logger.debug(f"Input net charge = {net_charge}")
    rounded_charge = int(round(net_charge))
    logger.debug(f"Rounding net charge to {rounded_charge}")

    return rounded_charge