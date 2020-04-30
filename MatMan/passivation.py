"""
Here all the passivation in volved things can be found
"""

def passivate_zinc_blende_slab(slab, passivant):
    """
    Passivates zincblende slabs (only conventional unit cells) with the type of atoms given in the inputs
    Inputs:
        slab: Atoms object which is a conventional unit cell of the slab that you would like to passivate
        passivant: String object with the species that you would like to passivate the slab with
    return: Atoms object
    """
    slab = slab.copy()  # To prevent any previous instances lurking
    slab *= (1, 1, 2)  # we are here doubling the slab to get thorse top atoms that we can convert to other atoms.
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

    return slab