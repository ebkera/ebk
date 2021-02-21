"""
MatMan stands for Material Manager. 
It uses the ASE library to manage some of the material that we will use
"""

def make_inversion_symmetric(atoms):
    """
    This will make the cell inversion symmetric. 
    Funtionality till now includes only inversion symmetry along the z axis.
        What this means is that all Z>0 atoms will be rewritten so that we have inversion symmetry which again means that
        the surface of interest should be oriented in the z axis
    """
    # print(atoms.cell)
    z_axis_length = atoms.cell[2][2]
    atoms.cell[2][2] = 2*z_axis_length

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

    # for x in range(len(self.atoms), -1, -1):
    #     if x in zg0:
    #         if x < self.sb: self.sb = self.sb - 1
    #         elif x == self.sb: self.sb = "Has been deleted"
    #         del self.atoms[x]
