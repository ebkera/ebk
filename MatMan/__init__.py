"""
MatMan stands for Material Manager. 
It uses the ASE library to manage some of the material that we will use
"""

def make_inversion_symmetric(atoms, x_offset = 0, y_offset = 0):
    """
    This will make the cell inversion symmetric. 
    Funtionality till now includes only inversion symmetry along the z axis.
        What this means is that all Z>0 atoms will be rewritten so that we have inversion symmetry which again means that
        the surface of interest should be oriented in the z axis
    """
    # print(atoms.cell)
    z_axis_length = atoms.cell[2][2]
    # atoms.cell[2][2] = 2*z_axis_length
    # atoms.cell[1][1] += x_offset
    # atoms.cell[2][2] += y_offset

    # Lets move all atoms downwards by the amount of z axis length of the cell
    for i,v in enumerate(atoms):
        v.position[2] -= z_axis_length

    zg0 = []  # Atoms that are at z>0
    ze0 = []  # Atoms that are at z=0
    zl0 = []  # Atoms that are at z<0

    for i,v in enumerate(atoms):
        if v.position[2] > 0: zg0.append(i)
        elif v.position[2] == 0: ze0.append(i)
        elif v.position[2] < 0: zl0.append(i)

    # print(f"zg0: {zg0}")
    # print(f"ze0: {ze0}")
    # print(f"zl0: {zl0}")
    # print(len(self.atoms))

    for x in zl0:
        atoms.append(atoms[x])
        last_atom = len(atoms)-1
        for i in range(3):
            atoms[last_atom].position[i] = -atoms[last_atom].position[i]

    # This was when self had other reasons for keeping track of atoms
    # for x in range(len(self.atoms), -1, -1):
    #     if x in zg0:
    #         if x < self.sb: self.sb = self.sb - 1
    #         elif x == self.sb: self.sb = "Has been deleted"
    #         del self.atoms[x]

def find_cell(atoms, a0, **kwargs):
    """
    x, y, z will add extra vacuum
    LI.find_cell(atoms, a0, z=a0/4, x=0, y=0, xcell=a0/4, ycell=a0/4)
    """
    x = kwargs.get("x", 0)  # Any vaccuum that you want to insert
    y = kwargs.get("y", 0)
    z = kwargs.get("z", 0)
    x_offset = kwargs.get("x_offset", 0)
    y_offset = kwargs.get("y_offset", 0)

    max_x = 0
    max_y = 0
    max_z = 0
    min_x = 1000000
    min_y = 1000000
    min_z = 1000000
    for atom in atoms:
        if atom.position[0] > max_x: max_x = atom.position[0]
        if atom.position[1] > max_y: max_y = atom.position[1]
        if atom.position[2] > max_z: max_z = atom.position[2]
        if atom.position[0] < min_x: min_x = atom.position[0]
        if atom.position[1] < min_y: min_y = atom.position[1]
        if atom.position[2] < min_z: min_z = atom.position[2]

    # calculating the cell offsets in the x, y, z directions. Here we assume the ligand is centred at the origin to startwith
    # x_offset = a0/4-#atoms.atoms[atoms.sb].position[0]*2
    # y_offset = -a0/4-#atoms.atoms[atoms.sb].position[1]*2

    # These are the cell vectors
    x_cell = (a0 + x, 0, 0)
    y_cell = (0, a0 + y, 0)
    z_cell = (x_offset, y_offset, max_z-min_z + a0/4)  # Here we have to add a0/4 since there is no information on that.

    cell = [x_cell, y_cell, z_cell]
    translation_vec = [(x_cell[0]+x_offset)/2, (y_cell[1]+y_offset)/2, (z_cell[2])/2]
    translation_vec = [(x_cell[0]+x_offset)/2, (y_cell[1]+y_offset)/2, (z_cell[2])/2]
    for atom in atoms:
        atom.position = (atom.position[0] + translation_vec[0], atom.position[1] + translation_vec[1], atom.position[2] + translation_vec[2])
    atoms.set_cell(cell)
    return cell
