"""
Here all the passivation in volved things can be found
"""
from ase.build import cut

def passivate_zinc_blende_slab(slab, passivant):
    """
    Passivates zincblende slabs (only conventional unit cells) with the type of atoms given in the inputs
    Inputs:
        slab: Atoms object which is a conventional unit cell of the slab that you would like to passivate
        passivant: String object with the species that you would like to passivate the slab with
    return: Atoms object
    """
    slab = slab.copy()  # To prevent any previous instances lurking
    slab *= (1, 1, 2)  # we are here doubling the slab to get those top atoms that we can convert to other atoms.
    # slab.rotate(90, '-x')
    # print(f"Number of atoms before deletion: {len(slab)}")

    for x in range(15, 10, -1):  # Delteing the extra atoms of the species up on top of the initial slab
        del(slab[x])
    del(slab[9])
    # print(f"The slab contains : {len(slab)} atoms in total")

    #Next We will change the two new atoms to passivation atoms  type and introduce another two atoms that we will change 
    slab[8].symbol = f"{passivant}"
    slab[9].symbol = f"{passivant}"
    slab.append(f"{passivant}")
    slab.append(f"{passivant}")
    coords = slab.get_positions()
    # print(coords)
    n = 8
    d = 1
    # fist lets copy the old H atom coordinates to the newly added H atoms
    coords[n+2] = [coords[n][0], coords[n][1], coords[n][2]]
    coords[n+3] = [coords[n+1][0], coords[n+1][1], coords[n+1][2]]
    # Then lets change the coordinates so that they are changed by d angs down and to the sides
    for n in range(8,10):
        # Important, These below statements have to be in this order. If reversed this will change n and then will affect the subsequent statements
        coords[n+2] = [coords[n][0]+d, coords[n][1]-d, coords[n][2]-d]
        coords[n] = [coords[n][0]-d, coords[n][1]+d, coords[n][2]-d]
    slab.set_positions(coords)
    slab.edit()
    return slab


def passivate_zinc_blende_slab_2(slab, passivant, direction):
    """
    The idea is to have the slabs passivated with passivant in the direction.
    right now only works for the 110 direction
    default zinc blende sturcture should be input as slab if not in the 001 direction
    """
    import numpy as np
    
    if direction == (0,0,1):
        return passivate_zinc_blende_slab(slab, passivant)
    else:
        c = np.array(direction, dtype=float)
        b = np.array((0,0,1), dtype=float)
        a = np.cross(b, direction)
        y = cut(slab, a, b*(1/2), nlayers=5)
        y.edit()
        # coords = y.get_positions()
        
        d = 2  # displace by this amount
        xy = .75
        if passivant == "H":
            harcoded_atom_numbers = [18, 19, 14, 16]
            for num in harcoded_atom_numbers:
                y[num].symbol = f"{passivant}"
                # print("before")
                # print(y[num].position)
                # y[num].position = [y[num].position[0], y[num].position[1], y[num].position[2] - d]
                y[num].position = [y[num].position[0] - xy, y[num].position[1] - xy, y[num].position[2]]

                # # New replaced atoms since H has two coordinations to take care of now
                # # since we have displaced through the xy directions already as above we dont need to do the new atoms since we have already got the right xy
                # y.append(f"{passivant}")
                # y[-1].position = [y[num].position[0], y[num].position[1], y[num].position[2] + 2*d]
                # print("after")
                # print(y[num].position)
                # print(y[-1].position)
        
        
        elif passivant == "O":
            pass
            
        return y