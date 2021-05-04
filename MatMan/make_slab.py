"""
Uses ASE to make slabs:
    ASE does not support higher order surfaces by default but offers the ase.build.cut method. 
    We use this method to quickly make slabs of different surfaces orientations.
"""

# import ase.build.cut as cut
from ase.lattice.compounds import Zincblende  # For future use especially if required for HgTe
from ase.lattice.cubic import Diamond
import ase
import ebk.MatMan
from ebk.coordinateMaker import CoordinateMaker as i_s_a
from ebk import progress_bar
from ase import Atom
import numpy as np
# from ebk.MatMan.insert_ligands 

class MakeSlab():
    def __init__(self, *args, **kwargs):
        # print(kwargs)
        self.lattice_constant = self.a0 = kwargs["a0"]

    def load(self, atoms_object, *args, **kwargs):
        """Loads external structures from atoms like objects"""
        self.atoms = atoms_object

    def add_vacuum(self, vacuum):
        """ Adds the amount of vacuum in angstroms"""
        ase.build.add_vacuum(self.atoms, vacuum)

    def make_inversion_symmetric(self):
        ebk.MatMan.make_inversion_symmetric(self.atoms)
        self.center_to_cell()

    def center_to_cell(self):
        self.atoms.center()

    def passivate_zinc_blende_slab(self, bond_length, passivation_direction = ["x", "y", "z"]):
        """This method will add hydrogen atoms to surface atoms of the dot that has any dangling bonds"""
        self.identify_surface_atoms()
        self.passivation_bondlength = bond_length
        multi_factor = bond_length/(self.lattice_constant*0.4330127018922)  # Since the coordinates are in lattice constants where 0.4330127018922 is the length of a bond 

        # These are the sites of the nearest neighbours
        positive_set_pos = [[0.25,0.25,0.25],[-0.25,-0.25,0.25],[-0.25,0.25,-0.25],[0.25,-0.25,-0.25]]
        negative_set_pos = [[-0.25,-0.25,-0.25],[0.25,0.25,-0.25],[0.25,-0.25,0.25],[-0.25,0.25,0.25]]
        positive_set = [["+","+", "+"],["-","-", "+"],["-","+","-"],["+","-","-"]]
        negative_set = [["-", "-", "-"],["+", "+", "-"],["+", "-", "+"],["-", "+", "+"]]
        new_H_atoms = []

        hydrogenate_progress = progress_bar(len(self.surface_atoms_list))
        print(f"hydrogenate: Checking and passivating")
        extreme_coordinates = [[100000,0], [100000,0], [100000,0]]  # list of 3 lists (x,y,z): inner list is of len 2 with - and + extremums

        for i, v in enumerate(self.surface_atoms_list):
            # Here we are finding out the extreme atoms. This will help in determining what directions get passivated
            # hydrogenate_progress.get_progress(i)
            displacement_vectors = []  # This will store the displacement vectors between current atoms connected to "atom" which is the site where we want to hydrogenate
            # This will store all the new Hydrogen atoms that will be introduced and we can add all of them at once at the end
            # Finding the current neighbours for the site and also the displacement vectors for each site
            atom = self.atoms[v].position
            for direction in [0,1,2]:
                # for sign in [0,1]: # where 0 and 1 are for - and +
                if atom[direction]<extreme_coordinates[direction][0]: extreme_coordinates[direction][0] = atom[direction]
                if atom[direction]>extreme_coordinates[direction][1]: extreme_coordinates[direction][1] = atom[direction]

        for i, v in enumerate(self.surface_atoms_list):
            # print(self.surface_atoms_list)
            hydrogenate_progress.get_progress(i)
            displacement_vectors = []  # This will store the displacement vectors between current atoms connected to "atom" which is the site where we want to hydrogenate
            # This will store all the new Hydrogen atoms that will be introduced and we can add all of them at once at the end
            # Finding the current neighbours for the site and also the displacement vectors for each site
            atom = self.atoms[v].position

            # Cheking for passivation directions
            continue_flag = False
            for direc in passivation_direction:
                if   direc.lower() == "x":
                    if atom[0] in extreme_coordinates[0]: pass
                    else: continue_flag = True
                elif direc.lower() == "y":
                    if atom[1] in extreme_coordinates[1]: pass
                    else: continue_flag = True
                elif direc.lower() == "z":
                    if atom[2] in extreme_coordinates[2]: pass
                    else: continue_flag = True
                else:
                    continue_flag = True

            # This will ensure that only the required directions are passivated
            if continue_flag: continue

            # for direction in [0,1,2]:
            #     # for sign in [0,1]: # where 0 and 1 are for - and +
            #     if atom[direction]<extreme_coordinates[direction][0]: extreme_coordinates[direction][0] = atom[direction]
            #     if atom[direction]>extreme_coordinates[direction][1]: extreme_coordinates[direction][1] = atom[direction]

            # print(atom)
            # print(extreme_coordinates)
            for atom_i in self.atoms:
                atom2 = atom_i.position
                displacement_vector_x = atom2[0] - atom[0]
                displacement_vector_y = atom2[1] - atom[1]
                displacement_vector_z = atom2[2] - atom[2]
                displacement_vector_mag2 = abs(displacement_vector_x) ** 2 + abs(displacement_vector_y) ** 2 + abs(displacement_vector_z) ** 2
                # print(atom, atom2)
                # print(displacement_vector_mag2, 0.1875 * self.a0**2, displacement_vector_mag2 - 0.1875 * self.a0 **2)
                if abs(displacement_vector_mag2 - 0.1875 * self.a0**2) <= 0.00000000001:
                    # print("found nn")
                    if float(displacement_vector_x) < 0: displacement_vector_x = "-"
                    else: displacement_vector_x = "+"
                    if float(displacement_vector_y) < 0: displacement_vector_y = "-"
                    else: displacement_vector_y = "+"
                    if float(displacement_vector_z) < 0: displacement_vector_z = "-"
                    else: displacement_vector_z = "+"
                    displacement_vector = [displacement_vector_x, displacement_vector_y, displacement_vector_z]
                    displacement_vectors.append(displacement_vector)
            current_missing_neighbours = []
            type_of_neighbours = ""
            # print(f"displacement vectors {displacement_vectors}")
            for i, v in enumerate(positive_set):
                type_of_neighbours = "positive"
                if v not in displacement_vectors:
                    current_missing_neighbours.append(i)
            if len(current_missing_neighbours) == 4: # if the count here is 4 then obviously we have to use the negative set
                current_missing_neighbours = []
                type_of_neighbours = "negative"
                for i, v in enumerate(negative_set):
                    if v not in displacement_vectors:
                        current_missing_neighbours.append(i)
            if type_of_neighbours == "positive": passivation_atoms_to_add_pos = positive_set_pos
            else: passivation_atoms_to_add_pos = negative_set_pos
            positions=[0,0,0]
            # print(current_missing_neighbours, passivation_atoms_to_add_pos)
            # print(type_of_neighbours)
            # atom

            for x in current_missing_neighbours:
                # Now we prevent from making towards the wrong direction 
                continue_flag = False
                # print(atom, extreme_coordinates, passivation_atoms_to_add_pos[x])
                for direc in passivation_direction:
                    if   direc.lower() == "x":
                        if abs(atom[0] - extreme_coordinates[0][0]) < 0.00001 and passivation_atoms_to_add_pos[x][0] < 0: pass
                        elif abs(atom[0] - extreme_coordinates[0][1]) < 0.00001 and passivation_atoms_to_add_pos[x][0] > 0: pass
                        else: continue_flag = True
                    elif direc.lower() == "y":
                        if abs(atom[1] - extreme_coordinates[1][0]) < 0.00001 and passivation_atoms_to_add_pos[x][1] < 0: pass
                        elif abs(atom[1] - extreme_coordinates[1][1]) < 0.00001 and passivation_atoms_to_add_pos[x][1] > 0: pass
                        else: continue_flag = True
                    elif direc.lower() == "z":
                        if abs(atom[2] - extreme_coordinates[2][0]) < 0.00001 and passivation_atoms_to_add_pos[x][2] < 0: pass
                        elif abs(atom[2] - extreme_coordinates[2][1]) < 0.00001 and passivation_atoms_to_add_pos[x][2] > 0: pass
                        else: continue_flag = True
                    else:
                        continue_flag = True

                # print(continue_flag)

                if continue_flag: continue

                positions[0] = passivation_atoms_to_add_pos[x][0]*bond_length*(np.cos(np.pi/4))**2/0.25 + atom[0]
                positions[1] = passivation_atoms_to_add_pos[x][1]*bond_length*(np.cos(np.pi/4))**2/0.25 + atom[1]
                positions[2] = passivation_atoms_to_add_pos[x][2]*bond_length*(np.cos(np.pi/4))**2/0.25 + atom[2]
                new_H_atoms.append(Atom('H', position=positions))
        # self.atoms.append(atom_to_append)

        # # self.no_of_passivation_atoms += len(new_H_atoms)
        for x in new_H_atoms:
            self.atoms.append(x)
        print("hydrogenate: Successfully passivated with hydrogen")

    def identify_surface_atoms(self):
        """Here we are just recognizing the sites that are not passivatied and have dangling bonds after trimming"""
        print(f"identify_surface_atoms: Started")
        surface_counter = 0
        length = len(self.atoms)
        # surface_progress = progress_bar(length)
        self.surface_atoms_list = []
        for i, atom in enumerate(self.atoms):
            # surface_progress.get_progress(i)
            nearest_neigbours = 0
            for atom2 in self.atoms:
                if abs(abs(atom.position[0] - atom2.position[0]) - 0.25*self.a0) <= 0.00000000001:
                    # print(f"{atom.position[0] - atom2.position[0]} vs {0.25*self.a0}" )
                    if abs(abs(atom.position[1] - atom2.position[1]) - 0.25*self.a0) <= 0.00000000001:
                        if abs(abs(atom.position[2] - atom2.position[2]) - 0.25*self.a0) <= 0.00000000001:
                            nearest_neigbours +=1
                if nearest_neigbours == 4:
                    break

            if nearest_neigbours != 4:
                surface_counter += 1
                self.surface_atoms_list.append(i)

        print(f"identify_surface_atoms: Done: Checked: {i+1} atoms           ")  # White space to prevent ghosting
        # print(self.surface_atoms_list)
        return surface_counter

class Diamond_bulk():
    def __init__(self, a0, nlayers, *args, **kwargs):
        self.atoms = Diamond(symbol="Sn", latticeconstant=a0, pbc=(1,1,1))
        super().__init__()

# class Diamond100(MakeSlab):
#     """This works but uses the cut method"""
#     def __init__(self, a0, nlayers, *args, **kwargs):
#         self.atoms = Diamond(symbol="Sn", latticeconstant=a0, pbc=(1,1,1))
#         self.atoms = ase.build.cut(self.atoms, a=(1, 0, 0), b=(0, 1, 0), c=None, clength=None, origo=(0, 0, 0), nlayers=nlayers, extend=1.0, tolerance=0.01, maxatoms=None)
#         super().__init__()

class Diamond100(MakeSlab):
    def __init__(self, a0, nlayers, *args, **kwargs):
        self.atoms = Diamond(symbol="Sn", latticeconstant=a0, directions=[[1,0,0], [0,1,0], [0,0,1]], size=(1,1,nlayers), pbc=(1,1,0), miller=[None, None, [0, 0, 1]])
        super().__init__(**{"a0":a0})

class Diamond110(MakeSlab):
    """Still under construction"""
    def __init__(self, a0, nlayers, *args, **kwargs):
        self.atoms = Diamond(symbol="Sn", latticeconstant=a0, directions=[[1,-1,0], [0,0,-1], [1,1,0]], size=(1,1,nlayers), pbc=(1,1,0), miller=[None, None, [1,1,0]])
        # self.atoms = ase.build.cut(self.atoms, a=(1, 0, 0), b=(0, 1, 0), c=None, clength=None, origo=(0, 0, 0), nlayers=nlayers, extend=1.0, tolerance=0.01, maxatoms=None)
        super().__init__(**{"a0":a0})

class Diamond111(MakeSlab):
    """Creates Slabs in with the 111 direction now lying along the z axis"""
    def __init__(self, a0, nlayers, *args, **kwargs):
        self.atoms = Diamond(symbol="Sn", latticeconstant=a0, directions=[[1,-1,0], [1,1,-2], [1,1,1]], size=(1,1,nlayers), pbc=(1,1,0), miller=[None, None, [1,1,1]])
        super().__init__(**{"a0":a0})

class Diamond210(MakeSlab):
    """Creates Slabs in with the 210 direction now lying along the z axis"""
    def __init__(self, a0, nlayers, *args, **kwargs):
        self.atoms = Diamond(symbol="Sn", latticeconstant=a0, directions=[[1,-2,0], [0,0,-1], [2,1,0]], size=(1,1,nlayers), pbc=(1,1,0), miller=[None, None, [2,1,0]])
        super().__init__(**{"a0":a0})
        
class Diamond211(MakeSlab):
    """Creates Slabs in with the 211 direction now lying along the z axis"""
    def __init__(self, a0, nlayers, *args, **kwargs):
        self.atoms = Diamond(symbol="Sn", latticeconstant=a0, directions=[[1,-2,0], [0,0,-1], [2,1,0]], size=(1,1,nlayers), pbc=(1,1,0), miller=[None, None, [2,1,1]])
        super().__init__(**{"a0":a0})