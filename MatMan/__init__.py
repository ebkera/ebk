"""
MatMan stands for Material Manager. 
It uses the ASE library to manage some of the material that we will use
"""

def make_inversion_symmetric(atoms, duplicate_z = "-"):
    """
    This will make the cell inversion symmetric. 
    inputs:
        duplicate z = (- or +) THis will make make either the negative or positive part of the cell as the original before inversion
    Funtionality till now includes only inversion symmetry along the z axis.
        What this means is that all Z>0 atoms will be rewritten so that we have inversion symmetry which again means that
        the surface of interest should be oriented in the z axis
    """
    zg0 = []  # Atoms that are at z>0
    ze0 = []  # Atoms that are at z=0
    zl0 = []  # Atoms that are at z<0

    for i,v in enumerate(atoms):
        if v.position[2] > 0: zg0.append(i)
        elif v.position[2] == 0: ze0.append(i)
        elif v.position[2] < 0: zl0.append(i)

    if duplicate_z == "+":
        atoms_to_invert = zg0
        atoms_to_delete = zl0
    else:
        atoms_to_invert = zl0
        atoms_to_delete = zg0

    for x in atoms_to_invert:
        atoms.append(atoms[x])
        last_atom = len(atoms)-1
        for i in range(3):
            atoms[last_atom].position[i] = -atoms[last_atom].position[i]

    for x in range(len(atoms), -1, -1):
        if x in atoms_to_delete:
            del atoms[x]

def find_cell(atoms, a0, slab_atoms_species="Sn", **kwargs):
    """
    This function as of now is intended for use with setups that contain slabs and ligands.
    It looks at the off set in the slab atoms on the either side of the ligand and looks to compensate for the change by changing the c axis
    """

    # We will first try to find the extreme slab atoms by looking at slab atoms extreme coordinates
    # Have not used the get_extreme_coordinates method since it might pick up ligand atoms (non slab atoms)
    max_z = 0
    min_z = 1000000 
    x_mf = kwargs.get("x_multiplication_factor", 1)
    y_mf = kwargs.get("y_multiplication_factor", 1)
    z_mf = kwargs.get("z_multiplication_factor", 1)

    for atom in atoms:
        if atom.symbol == slab_atoms_species:
            if atom.position[2] > max_z: max_z = atom.position[2]
            if atom.position[2] < min_z: min_z = atom.position[2]
        # print(atom.symbol)

    for atom in atoms:
        if atom.position[2] == max_z and atom.symbol == slab_atoms_species:
            max_z_atom_index = atom.index
        if atom.position[2] == min_z and atom.symbol == slab_atoms_species:
            min_z_atom_index = atom.index

    # z_length = extreme_coordinates[2][1] - extreme_coordinates[2][0]
    z_length = max_z - min_z 
    x_spill_over = atoms[max_z_atom_index].position[0] - atoms[min_z_atom_index].position[0] - a0/4
    y_spill_over = atoms[max_z_atom_index].position[1] - atoms[min_z_atom_index].position[1] + a0/4

    # These are the cell vectors
    # Original for slabs going as 1 UC 
    x_cell = (a0*x_mf, 0, 0)
    y_cell = (0, a0*y_mf, 0)

    # x_cell = (a0*x_mf, a0*y_mf, 0)
    # y_cell = (0, a0*y_mf, 0)
    # x_cell = (0, a0*y_mf, a0*z_mf)
    # y_cell = (a0*x_mf,0, a0*z_mf)
    # x_cell = (0, a0*y_mf, 0)
    # y_cell = (a0*x_mf,0, 0)
    z_cell = (x_spill_over, y_spill_over, z_length + a0/4)  # Here we have to add a0/4 since there is no information on that.
    # z_cell = (x_spill_over*x_mf, y_spill_over*y_mf, z_length + a0/4)  # Here we have to add a0/4 since there is no information on that.

    cell = [x_cell, y_cell, z_cell]
    atoms.set_cell(cell)
    atoms.center()
    return cell

def get_extreme_coordinates(atoms: 'atomsobject') -> list:
    """This function outputs extreme corrdinates of a atoms like object"""
    extreme_coordinates = [[100000,0], [100000,0], [100000,0]]  # list of 3 lists (x,y,z): inner list is of len 2 with - and + extremums
    for i, v in enumerate(atoms):
        atom = v.position
        for direction in [0,1,2]:
            # for sign in [0,1]: # where 0 and 1 are for - and +
            if atom[direction]<extreme_coordinates[direction][0]: extreme_coordinates[direction][0] = atom[direction]
            if atom[direction]>extreme_coordinates[direction][1]: extreme_coordinates[direction][1] = atom[direction]
    return extreme_coordinates

def get_center(atoms: 'atomsobject') -> list:
    """This centre of the atoms (not the cell)"""
    extreme_coordinates = get_extreme_coordinates(atoms)
    return([extreme_coordinates[0][1]-extreme_coordinates[0][0], extreme_coordinates[1][1]-extreme_coordinates[1][0], extreme_coordinates[2][1] - extreme_coordinates[2][0]])

def make_common_centre(atoms1: 'atomsobject', atoms2: 'atomsobject') -> list:
    """atoms2 will be shifted so that centre is at atoms1"""
    
    atom1_extreme_coordinates = get_extreme_coordinates(atoms1)
    atoms2_extreme_coordinates = get_extreme_coordinates(atoms2)

    centre_line_atom1 = []
    centre_line_atom2  = []
    offset = []

    for i in range(0,3):
        centre_line_atom1.append((atom1_extreme_coordinates[i][0] + atom1_extreme_coordinates[i][1])/2)
        centre_line_atom2.append((atoms2_extreme_coordinates[i][0] + atoms2_extreme_coordinates[i][1])/2)
        offset.append(centre_line_atom1[i] - centre_line_atom2[i])

    for i,v in enumerate(atoms2):
        for j in range(0,3):
            v.position[j] += offset[j]

    return atoms2

def is_inversion_symmetric(atoms) -> bool:
    """

    THis is under construction and is not in wrorking condition


    This function will check if the atoms type object is inversion symmetric.
    The function checks to see if the atoms object is actually inversion symmetric wrt to the cell
    """
    extreme_coordinates = get_extreme_coordinates(atoms)
    inverted_positions = []
    center = get_center(atoms)
    cell = atoms.get_cell()
    print(cell[2][2])

    for i,v in enumerate(atoms):
        inverted_positions.append((-v.position[0], -v.position[1], -v.position[2]))
        for j,w in enumerate(atoms):
            if ((w.position[0] == inverted_positions[j][0]) and (w.position[1] == inverted_positions[1]) and (w.position[2] == inverted_positions[2])):
                return False
    return True

