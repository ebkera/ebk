"""
Here all the passivation in volved things can be found
"""
from ase.build import cut
import numpy as np

def passivate_zincblende_slab_001(slab, passivant):
    """
    Passivates zincblende slabs (Only conventional unit cell distances accepted in the non 001 direction) with the type of atoms given in the inputs
    Inputs:
        slab: Atoms object which is a conventional unit cell of the slab that you would like to passivate
              Can be any number of unit cells in teh z direction.
        passivant: String object with the species that you would like to passivate the slab with
    return: Atoms object
    """
    slab = slab.copy()  # To prevent any previous instances lurking
    # we need to get atoms that match these ones
    move_up = [0,2]
    # move_down = [5,7]
    move_down = [len(slab.get_positions())-1, len(slab.get_positions())-3]
    xy = 1
    d = slab.get_cell()[2][2]
    for v in move_up:
        if passivant == "H":
            slab.append(f"{passivant}")
            slab[-1].position = np.array([slab[v].position[0] + xy, slab[v].position[1] - xy, slab[v].position[2] + d], dtype=float)
            slab.append(f"{passivant}")
            slab[-1].position = np.array([slab[v].position[0] - xy, slab[v].position[1] + xy, slab[v].position[2] + d], dtype=float)
    for v in move_down:
        if passivant == "H":
            slab.append(f"{passivant}")
            slab[-1].position = np.array([slab[v].position[0] + xy, slab[v].position[1] - xy, slab[v].position[2] - d], dtype=float)
            slab.append(f"{passivant}")
            slab[-1].position = np.array([slab[v].position[0] - xy, slab[v].position[1] + xy, slab[v].position[2] - d], dtype=float)
    # slab.edit()
    return slab

def passivate_zinc_blende_slab_nonprimitive(slab, passivant, direction):
    """
    The idea is to have the slabs passivated with passivant in the direction.
    right now only works for the 110 direction
    default zinc blende sturcture should be input as slab if not in the 001 direction
    """
    import numpy as np
    
    if direction == (0,0,1):
        return passivate_zincblende_slab_001(slab, passivant)
    else:
        c = np.array(direction, dtype=float)
        b = np.array((0,0,1), dtype=float)
        a = np.cross(b, direction)
        y = cut(slab, a, b, nlayers=6)
        # y.edit()
        # coords = y.get_positions()
        
        # Fixing the amount that you want to displace by
        d = 2  # displace by this amount
        xy = .75

        if passivant == "H":
            atoms_up = [20, 15, 23, 22]  # These are atoms that were coreated using an extra layer at the top
            atoms_down = [0, 6, 3, 2]  # These are atoms that were coreated using an extra layer at the bottom
            for num in atoms_up:
                y[num].symbol = f"{passivant}"
                # Since in the 110 direction we have only the one bond to take care of in zincblende we need only replace the atoms with all of these
                y[num].position = [y[num].position[0], y[num].position[1], y[num].position[2]]
            for num in atoms_down:
                y[num].symbol = f"{passivant}"
                # Since in the 110 direction we have only the one bond to take care of in zincblende we need only replace the atoms with all of these
                y[num].position = [y[num].position[0], y[num].position[1], y[num].position[2]]
        elif passivant == "O":
            pass
        return y

def passivate_zinc_blende_slab(slab, passivant, direction, n_bilayers, primtive = True):
    """
    The idea is to have the slabs passivated with passivant in the direction.
    right now only works for the 110 direction
    default zinc blende sturcture should be input as slab if not in the 001 direction
    n_bilayers is the number of bilayers in the slab
    """
    import numpy as np

    if direction == (0,0,1):
        return passivate_zincblende_slab_001(slab, passivant)
    else:
        c = np.array(direction, dtype=float)
        b = np.array((0,0,1), dtype=float)
        a = np.cross(b, direction)

        # here we set the scale of the xy directions
        if primtive == True: xy_scale = 1
        if primtive == False: xy_scale = 2

        y = cut(slab, xy_scale*a/2, b, nlayers=2*n_bilayers+2)
        # y.edit()
        # coords = y.get_positions()
        
        # Fixing the amount that you want to displace by
        d = 2  # displace by this amount
        xy = .75

        if passivant == "H":
            if primtive == True:
                atoms_up = [-1, -2]  # These are atoms that were coreated using an extra layer at the top
                atoms_down = [0, 1]  # These are atoms that were coreated using an extra layer at the bottom
            if primtive == False:
                atoms_up = [-4, -9, -1, -2]  # These are atoms that were coreated using an extra layer at the top
                atoms_down = [0, 6, 3, 2]  # These are atoms that were coreated using an extra layer at the bottom                
            for num in atoms_up:
                y[num].symbol = f"{passivant}"
                # Since in the 110 direction we have only the one bond to take care of in zincblende we need only replace the atoms with all of these
                y[num].position = [y[num].position[0], y[num].position[1], y[num].position[2]]
            for num in atoms_down:
                y[num].symbol = f"{passivant}"
                # Since in the 110 direction we have only the one bond to take care of in zincblende we need only replace the atoms with all of these
                y[num].position = [y[num].position[0], y[num].position[1], y[num].position[2]]                
        elif passivant == "O":
            pass
        return y