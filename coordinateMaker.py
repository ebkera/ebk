"""
This code has the following functionality
Can make diamond/zinc blende material to required size blocks
Can trim these blocks to make sperical quantum dots
Can passivate these dots with Hydrogen atoms (can change H bond length)
Can output coordinates in multimple formats including .xyz and .fdf formats
Can output file with nearest neighbours and next nearest neighbours for use in any Tight Binding code
"""
import os
import sys
import subprocess
import time
import matplotlib

matplotlib.use('Agg')  # no UI backend required if working in the wsl without a UI
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

def radians_to_degrees(rad):
    """Converts radians into degrees"""
    return rad*180/np.pi

def triatom_angle(atom, atom1, atom2):
    """Calculates that angle that 'atom' makes with 'atom1' and 'atom2'"""
    v1 = [atom1[0]-atom[0],atom1[1]-atom[1],atom1[2]-atom[2]]
    v2 = [atom2[0]-atom[0],atom2[1]-atom[1],atom2[2]-atom[2]]
    return radians_to_degrees(np.arccos((v1[0]*v2[0]+v1[1]*v2[1]+v1[1]*v2[1])/(np.linalg.norm(v1)*np.linalg.norm(v2))))

def invert_coordinates(input_array):
    """Does inversion symmetry on all the coordinates of the given NP"""
    output_array = []
    for x in input_array:
        output_array.append([-x[0], -x[1], -x[2], x[3]])
    return output_array

def cartesian_to_polar(u2,u1):
    """This fuction converts the vector that points from u1 to u2 from cartesian to spherical polar"""
    x = u2[0]-u1[0]
    y = u2[1]-u1[1]
    z = u2[2]-u1[2]
    # Calculating r
    r = np.sqrt(x**2+y**2+z**2)
    # Calculating theta
    if z == 0.0:
        theta = 90
    else:
        theta = radians_to_degrees( np.arctan(np.sqrt((x**2+y**2)/z**2)))
    # Calculating phi
    if x == 0.0 and y > 0.0:
        phi = 90.0
    elif x == 0.0 and y < 0.0:
        phi = 90.0
    elif x < 0 and y > 0.0 or x < 0 and y < 0.0:
            phi = radians_to_degrees(np.arctan(y/x) + np.pi)
    else:
        phi = radians_to_degrees(np.arctan(y/x))
    return([r,theta,phi])

class DotDiamond():
    def __init__(self, radius, coordinate_format, lattice_constant, replicate):
        """ The radius is the half-length if a side of a cube that will be created.
        The cutoff is the actual radius of the dot we use the cutoff to trim the dots that are further than the cutoff length"""
        self.conventional_cell = []
        # This is the conventional (FCC) unit cell for a diamond lattice
        self.conventional_cell.append([0.000, 0.000, 0.000, "A"])  # Index:0 atom at 000 (fcc_a)
        self.conventional_cell.append([0.250, 0.250, 0.250, "C"])  # Index:1 basis atom connected to (fcc_b)
        # self.conventional_cell.append([0.000,0.500,0.500]) # Index:2 three atoms forming the lattice vectors (fcc_a). # This atom can be removed since slotting another box will duplicate this atom
        # self.conventional_cell.append([0.500,0.000,0.500]) # Index:3 # This atom can be removed since slotting another box will duplicate this atom
        # self.conventional_cell.append([0.500,0.500,0.000]) # Index:4 # This atom can be removed since slotting another box will duplicate this atom
        # self.conventional_cell.append([1.000,1.000,0.000]) # Index:4,5 last atom on the xy plane (fcc_a)
        self.conventional_cell.append([0.750, 0.750, 0.250, "C"])  # Index:5,6 last atom on the xy_.25 plane (fcc_b)
        self.conventional_cell.append([0.500, 1.000, 0.500, "A"])  # Index:6,7 last two atoms on the 0.5 xy plane (fcc_a)
        self.conventional_cell.append([1.000, 0.500, 0.500, "A"])  # Index:8
        self.conventional_cell.append([0.250, 0.750, 0.750, "C"])  # Index:7,9 The last two atoms on the 0.75 xy plane (fcc_b)
        self.conventional_cell.append([0.750, 0.250, 0.750, "C"])  # Index:8,10
        # self.conventional_cell.append([0.000,1.000,1.000]) # Index:9,11 three atoms forming z=1 xy plane (fcc_a) # This atom can be removed since slotting another box will duplicate this atom
        # self.conventional_cell.append([1.000,0.000,1.000]) # Index:10,12 # This atom can be removed since slotting another box will duplicate this atom
        self.conventional_cell.append([0.500, 0.500, 1.000, "A"])  # Index:11,13

        self.radius = radius
        self.replicate = replicate
        self.lattice_constant = lattice_constant
        self.coordinate_format = coordinate_format
        self.start_at = [0.000, 0.000, 0.000]
        self.now_at = self.start_at.copy()
        self.basis = [0.250, 0.250, 0.250]
        self.lattice_vectors = [[0.000, 0.500, 0.500], [0.500, 0.000, 0.500], [0.500, 0.500, 0.000]]
        self.Number_of_atoms_in_cell = 8
        self.super_cell = []  # This is the supercell before trimming
        self.anion_count = 0  # Keeping track of the anion and cations in the dot
        self.cation_count = 0  # Keeping track of the anion and cations in the dot
        self.evenize = False  # Keeping track of whether the dot was evenized to have equal number of anions and cations
        self.finalcell = []  # This is the final dot. All passivation and evenizing will be included

        # Building the super_cell
        # We start with the first octant (+,+,+)
        # First we make a row in the x direction that has the requested radius(number of boxes)
        # We always iterate upto radius + 1 so that we don't miss any atoms that might have radii that are at the edge
        # We start the row with 0 since we want to start at zero
        # Then for the rest we have to start at 1 since we don't want the row we just did to be counted or the xy square we just did to be counted
        for step in range(0, self.radius + 1):
            for coordinate in self.conventional_cell:
                new_vector_x = [coordinate[0] + float(step), coordinate[1], coordinate[2], coordinate[3]]
                self.super_cell.append(new_vector_x)
        # Then we replicated those rows in the y direction
        full_x_row = self.super_cell.copy()  # Saving the row of the super cell
        for step in range(1, self.radius + 1):
            for coordinate in full_x_row:
                new_vector_y = [coordinate[0], coordinate[1] + float(step), coordinate[2], coordinate[3]]
                self.super_cell.append(new_vector_y)
        # Then we replicate this base into stacks in the z direction
        full_xy_square = self.super_cell.copy()  # Saving the xy base in the first quadrant of the super cell
        for step in range(1, self.radius + 1):
            for coordinate in full_xy_square:
                new_vector_z = [coordinate[0], coordinate[1], coordinate[2] + float(step), coordinate[3]]
                self.super_cell.append(new_vector_z)

        if self.replicate == True:
            # Then we replicate this into the (-,+,+) octant
            full_ppp_cube = self.super_cell.copy()  # Saving the +++ cube
            for coordinate in full_ppp_cube:
                new_vector = [coordinate[0] - float(self.radius + 1), coordinate[1], coordinate[2], coordinate[3]]
                self.super_cell.append(new_vector)
            # Then we replicate this into the (Z>0) volume
            full_0pp_cube = self.super_cell.copy()  # Saving the +++ cube
            for coordinate in full_0pp_cube:
                new_vector = [coordinate[0], coordinate[1] - float(self.radius + 1), coordinate[2], coordinate[3]]
                self.super_cell.append(new_vector)
            # Then we replicate this into the full volume
            full_zgt0_cube = self.super_cell.copy()  # Saving the +++ cube
            for coordinate in full_zgt0_cube:
                new_vector = [coordinate[0], coordinate[1], coordinate[2] - float(self.radius + 1), coordinate[3]]
                self.super_cell.append(new_vector)  # Then we replicated those rows in the y direction
        self.finalcell = self.super_cell

    def T2SL(self, monolayers):
        """Expected that the input be a non replicated slab. All coordinates are postive values"""
        self.monolayers = monolayers
        file = open("before", "w")
        file2 = open("after", "w")
        file3 = open("in", "w")
        file4 = open("out", "w")
        for cell in self.finalcell:
            file3.write(f"{cell}\n")

        for cell in self.finalcell:
            for x in range(0,3):
                file.write(f"{cell}\n")
                if cell[x] >= 2*0.25*self.monolayers[x]:
                    # print(f"Has to be removed since {cell[x]} >= {self.monolayers[x]}")
                    file2.write(f"{cell}\n")
                    try:
                        # This try is because if the cell was already removed due to another coordinate being too large they it will raise an error
                        self.finalcell.remove(cell)
                        print(f'Cell {cell} made it')
                    except:
                        print(f'Cell {cell} made it to except')
                else:
                    print(f'Cell {cell} did not make it')
                    
        for cell in self.finalcell:
            file4.write(f"{cell}\n")
        
        file.close()
        file2.close()
        file3.close()
        file4.close()

    def trim_to_dot(self, cut_off, evenize): 
        """This method trims the initial cube into a ball"""
        self.finalcell = []  # This is since we have put self.finalcell = self.super_cell if not trimmed
        self.cut_off = cut_off
        # Removing atoms if position is over the radius
        surface_counter = 0
        print(f"trim_to_dot: Number of atoms before trimming: {len(self.super_cell)}")
        for atom in self.super_cell:
            displacement2 = atom[0] ** 2 + atom[1] ** 2 + atom[2] ** 2
            if displacement2 < (self.cut_off) ** 2:  # This is where any atom that is inside the cutoff is saved
                self.finalcell.append(atom)

        # Here we are just recognizing the sites that are not passivatied and have dangling bonds after trimming
        for atom in self.finalcell:
            nearest_neigbours = 0
            for atom2 in self.finalcell:
                vec = abs(atom[0] - atom2[0]) ** 2 + abs(atom[1] - atom2[1]) ** 2 + abs(atom[2] - atom2[2]) ** 2
                if vec == 0.1875:
                    nearest_neigbours += 1

            if nearest_neigbours == 1 or nearest_neigbours == 2 or nearest_neigbours == 3:
                surface_counter += 1
                if atom[3] == "A":
                    atom[3] = "SA"
                elif atom[3] == "C":
                    atom[3] = "SC"

        # Count the number of anions and cations
        self.anion_count = 0
        self.cation_count = 0
        for x in self.finalcell:
            if x[3] == "A" or x[3] == "SA":
                self.anion_count += 1
            elif x[3] == "C" or x[3] == "SC":
                self.cation_count += 1

        # Here we do the evenization of the number of cations and the number of anions
        if evenize:
            self.evenize = True
            ac_difference = self.anion_count - self.cation_count
            print(f"trim_to_dot: The anion cation difference: {ac_difference}")
            i = 0
            if ac_difference > 0:
                while i < ac_difference:
                    for x in self.finalcell:
                        if x[3] == "SA":
                            self.finalcell.remove(x)
                            i+=1
                            break

            elif ac_difference < 0:
                while i < abs(ac_difference):
                    for x in self.finalcell:
                        if x[3] == "SC":
                            self.finalcell.remove(x)
                            i+=1
                            break

        # Now that we have removed surface atoms we need to find the new surface atoms so that we can passivate them as well
        # Here we are just recognizing the sites that are not passivatied and have dangling bonds after trimming
            surface_counter = 0   
            for atom in self.finalcell:
                nearest_neigbours = 0
                for atom2 in self.finalcell:
                    vec = abs(atom[0] - atom2[0]) ** 2 + abs(atom[1] - atom2[1]) ** 2 + abs(atom[2] - atom2[2]) ** 2
                    if vec == 0.1875:
                        nearest_neigbours += 1

                if nearest_neigbours == 1 or nearest_neigbours == 2 or nearest_neigbours == 3:
                    surface_counter += 1
                    if atom[3] == "A":
                        atom[3] = "SA"
                    elif atom[3] == "C":
                        atom[3] = "SC"

        # Count the number of anions and cations
            self.anion_count = 0
            self.cation_count = 0
            for x in self.finalcell:
                if x[3] == "A" or x[3] == "SA":
                    self.anion_count += 1
                elif x[3] == "C" or x[3] == "SC":
                    self.cation_count += 1


        self.no_of_surface_atoms = surface_counter
        self.no_of_atoms_after_trimming = len(self.finalcell)
        print(f"trim_to_dot: Number of atoms after trimming: {self.no_of_atoms_after_trimming}")

    def write_to_log(self):
        """Printing the final into an out file that contains the coordinates"""
        file_my = open("coordinates.log", "w+")
        for i in self.finalcell:
            for j in i:
                try:
                    file_my.write(str(j) + "  ")
                except:
                    file_my.write(j + "  ")
            file_my.write("\n")

        # Printing out general information
        # Creating a letex table
        file_my.write("\n\n##--Log of run in Latex format--(Copy and paste this into .tex file)--\n\n")
        file_my.write("\\begin{table}[h]\n")
        file_my.write("\\centering\n")
        file_my.write("\\begin{tabular}{ll}\n")
        file_my.write("Atom count for initial box    &: " + str(len(self.super_cell)) + " \\\\\n")
        file_my.write("Initial box length            &: " + str(self.radius*2.0*self.lattice_constant) + " (\\AA) (" + str(self.radius*2.0) + " conventional unit cells)\\\\\n")
        file_my.write("Coordinate Format             &: " + str(self.coordinate_format)+" \\\\\n")
        try:
            file_my.write("Atom count after trimming     &: " + str(self.no_of_atoms_after_trimming) + " \\\\\n")
            file_my.write("Cutoff parameter set to       &: " + str(self.cut_off) + " conventional unit cells\\\\\n")
            file_my.write("Anion/Cation Ratio            &: " + str(self.anion_count/self.cation_count) + "\\\\\n")
        except:
            file_my.write("% The initial box was never trimmed!\n")
        file_my.write("Anion Count                   &: " + str(self.anion_count) + "\\\\\n")
        file_my.write("Cation Count                  &: " + str(self.cation_count) + "\\\\\n")
        try:
            file_my.write("Surface sites                 &: " + str(self.no_of_surface_atoms) + "\\\\\n")
        except:
            file_my.write("% The initial box was never trimmed! - No Surface atoms!\n")
        try:
            file_my.write("Passivation atom count        &: " + str(self.no_of_atoms_for_passivation) + "\\\\\n")
            file_my.write("Passivation bond length       &: " + str(self.passivation_bondlength) + " (\\AA)\\\\\n")
        except:
            file_my.write("% The initial box was never pasivated! - No Passivation atoms!\n")
        try:
            file_my.write("Dot diameter                  &: " + str(self.cut_off*2*self.lattice_constant)+ " (\\AA) (" + str(self.cut_off*2)+ " conventional unit cells)\\\\\n")
            file_my.write("Dot diameter with passivation &: $\\approx$ " + str(self.cut_off*2*self.lattice_constant+self.passivation_bondlength)+ " (\\AA)\\\\\n")
        except:
            file_my.write("% The initial box was not trimmed! - No defined cutoff - No diameter or radius!\n")
        file_my.write("Total atoms                   &: " + str(len(self.finalcell)) + "\n")
        file_my.write("\\end{tabular}\n")
        file_my.write("\\caption{QD with " + str(len(self.finalcell)) + " total number of atoms} \\label{tab:QD" + str(len(self.finalcell)) + "}\n")
        file_my.write("\\end{table}\n")
        file_my.close()
        print("write_to_log: Successfully written to the log file")

    def number_atoms(self, type):
        """Numbers the atoms in the final cell list so that they know which atom is which"""
        counter = 1
        if type == "creation":
            for x in self.finalcell:
                x.append(counter)
                counter += 1
        print(f"number_atoms: Successfully numbered atoms in order of {type}")


    def hydrogenate(self, bond_length):
        """This method will add hydrogen atoms to surface atoms of the dot that has any dangling bonds"""
        self.passivation_bondlength = bond_length
        multi_factor = bond_length/(self.lattice_constant*0.4330127018922)  # Since the coordinates are in lattice constants where 0.4330127018922 is the length of a bond 
        # These are the sites of the nearest neighbours
        positive_set = [[0.25,0.25,0.25],[-0.25,-0.25,0.25],[-0.25,0.25,-0.25],[0.25,-0.25,-0.25]]
        negative_set = [[-0.25,-0.25,-0.25],[0.25,0.25,-0.25],[0.25,-0.25,0.25],[-0.25,0.25,0.25]]
        new_H_atoms = []
        for atom in self.finalcell:
            displacement_vectors = []  # This will store the displacement vectors between current atoms connected to "atom" which is the site where we want to hydrogenate
            # This will store all the new Hydrogen atoms that will be introduced and we can add all of them at once at the end
            # Taking only the surface atoms of the quantum dot
            if atom[3] == "SC" or atom[3] == "SA":
                # Finding the current neighbours for the site and also the displacement vectors for each site
                for atom2 in self.finalcell:
                    displacement_vector_x = atom2[0] - atom[0]
                    displacement_vector_y = atom2[1] - atom[1]
                    displacement_vector_z = atom2[2] - atom[2]
                    displacement_vector_mag2 = abs(displacement_vector_x) ** 2 + abs(displacement_vector_y) ** 2 + abs(displacement_vector_z) ** 2
                    if displacement_vector_mag2 == 0.1875:
                        displacement_vector = [displacement_vector_x,displacement_vector_y,displacement_vector_z]
                        displacement_vectors.append(displacement_vector)
                to_edit = []
                for v in positive_set:
                    if v not in displacement_vectors:
                        to_edit.append(v)
                if len(to_edit) == 4: # if the count here is 4 then obviously we have to use the negative set
                    to_edit = []
                    for v in negative_set:
                        if v not in displacement_vectors:
                            to_edit.append(v)
                for x in to_edit:
                    new_atom = [atom[0]+x[0]*multi_factor,atom[1]+x[1]*multi_factor,atom[2]+x[2]*multi_factor,"H"]
                    new_H_atoms.append(new_atom)
        self.no_of_atoms_for_passivation = len(new_H_atoms)
        for x in new_H_atoms:
            self.finalcell.append(x)
        print("hydrogenate: Successfully passivated with hydrogen")


    def oxygenate(self, bond_length):
        """This method will add Oxyen atoms to surface atoms of the dot that has dangling bonds. All dangling bonds will be passivated by one Oxygen atom in this method"""
        self.passivation_bondlength = bond_length
        multi_factor = bond_length/self.lattice_constant  # Since the coordinates are in lattice constants
        
        # These are the sites of the nearest neighbours
        positive_set = [[0.25,0.25,0.25],[-0.25,-0.25,0.25],[-0.25,0.25,-0.25],[0.25,-0.25,-0.25]]
        negative_set = [[-0.25,-0.25,-0.25],[0.25,0.25,-0.25],[0.25,-0.25,0.25],[-0.25,0.25,0.25]]
        new_O_atoms = []
        for atom in self.finalcell:
            displacement_vectors = []  # This will store the displacement vectors between current atoms connected to "atom" which is the site where we want to hydrogenate
            # This will store all the new Hydrogen atoms that will be introduced and we can add all of them at once at the end
            # Taking only the surface atoms of the quantum dot
            if atom[3] == "SC" or atom[3] == "SA":
                # Finding the current neighbours for the site and also the displacement vectors for each site
                for atom2 in self.finalcell:
                    displacement_vector_x = atom2[0] - atom[0]
                    displacement_vector_y = atom2[1] - atom[1]
                    displacement_vector_z = atom2[2] - atom[2]
                    displacement_vector_mag2 = abs(displacement_vector_x) ** 2 + abs(displacement_vector_y) ** 2 + abs(displacement_vector_z) ** 2
                    if displacement_vector_mag2 == 0.1875:
                        displacement_vector = [displacement_vector_x,displacement_vector_y,displacement_vector_z]
                        displacement_vectors.append(displacement_vector)
                to_edit = []
                for v in positive_set:
                    if v not in displacement_vectors:
                        to_edit.append(v)
                if len(to_edit) == 4: # if the count here is 4 then obviously we have to use the negative set
                    to_edit = []
                    for v in negative_set:
                        if v not in displacement_vectors:
                            to_edit.append(v)
                x_mag = 0
                y_mag = 0
                z_mag = 0
                for x in to_edit:  # Here we take the atoms that need to be passivated and then 
                    x_mag += x[0]
                    y_mag += x[1]
                    z_mag += x[2]
                to_edit = [[x_mag,y_mag,z_mag]]
                for x in to_edit:
                    new_atom = [atom[0]+x[0]*multi_factor,atom[1]+x[1]*multi_factor,atom[2]+x[2]*multi_factor,"H"]
                    new_O_atoms.append(new_atom)
        self.no_of_atoms_for_passivation = len(new_O_atoms)
        for x in new_O_atoms:
            self.finalcell.append(x)

    def write_to_fdf(self, zincblende, surface):
        # Printing the final into an out file that contains the coordinates
        file_fdf = open("coordinates.fdf", "w+")
        if self.coordinate_format == "Ang":
            file_fdf.write("AtomicCoordinatesFormat  Ang\n")
        elif self.coordinate_format == "Fractional":
            file_fdf.write("AtomicCoordinatesFormat  Fractional\n")
        file_fdf.write("%block AtomicCoordinatesAndAtomicSpecies\n")
        for i in self.finalcell:
            if self.coordinate_format == "Ang":
                file_fdf.write(
                    str(i[0] * self.lattice_constant) + "  " + str(i[1] * self.lattice_constant) + "  " + str(
                        i[2] * self.lattice_constant) + "  ")
            elif self.coordinate_format == "Fractional":
                file_fdf.write(str(i[0]) + "  " + str(i[1]) + "  " + str(i[2]) + "  ")
            # For writing the Type of atom
            if i[3] == "A":
                file_fdf.write("1")
            elif i[3] == "C" and zincblende == False:
                file_fdf.write("1")
            elif i[3] == "C" and zincblende == True:
                file_fdf.write("2")
            elif i[3] == "SA" or i[3] == "SC" and zincblende == False and surface == False:
                file_fdf.write("1")
            elif i[3] == "SA" or i[3] == "SC" and zincblende == False and surface == True:
                file_fdf.write("2")
            elif i[3] == "SA" or i[3] == "SC" and zincblende == True and surface == False:
                if i[3] == "SA":
                    file_fdf.write("1")
                elif i[3] == "SC":
                    file_fdf.write("2")
            elif i[3] == "SA" or i[3] == "SC" and zincblende == True and surface == True:
                file_fdf.write("3")
            elif i[3] == "H" and zincblende == False:
                file_fdf.write("2")
            elif i[3] == "H" and zincblende == True and surface == False:
                file_fdf.write("3")
            elif i[3] == "H" and zincblende == True and surface == True:
                file_fdf.write("4")
            file_fdf.write("\n")
        file_fdf.write("%endblock AtomicCoordinatesAndAtomicSpecies\n")
        print("write_to_fdf: Successfully written to fdf file")
        file_fdf.close()

    def write_to_fdf_zmatrix(self, zincblende, surface):
        # Printing the final into an out file that contains the coordinates in Z matrix format
        print(f"write_to_fdf_zmatrix: Make sure that the central atom is at the origin (cartesian: 0 0 0)")
        file_fdf = open("coordinates.fdf", "w+")
        file_fdf.write("%block Zmatrix\n")
        file_fdf.write("molecule\n")
        file_fdf.write("cartesian\n")
        for i in range(0,len(self.finalcell)):
        # for i in self.finalcell:
            # For writing the Type of atom
            if self.finalcell[i][3] == "A":
                file_fdf.write("1  ")
            elif self.finalcell[i][3] == "C" and zincblende == False:
                file_fdf.write("1  ")
            elif self.finalcell[i][3] == "C" and zincblende == True:
                file_fdf.write("2  ")
            elif self.finalcell[i][3] == "SA" or self.finalcell[i][3] == "SC" and zincblende == False and surface == False:
                file_fdf.write("1  ")
            elif self.finalcell[i][3] == "SA" or self.finalcell[i][3] == "SC" and zincblende == False and surface == True:
                file_fdf.write("2  ")
            elif self.finalcell[i][3] == "SA" or self.finalcell[i][3] == "SC" and zincblende == True and surface == False:
                if self.finalcell[i][3] == "SA":
                    file_fdf.write("1  ")
                elif self.finalcell[i][3] == "SC":
                    file_fdf.write("2  ")
            elif self.finalcell[i][3] == "SA" or self.finalcell[i][3] == "SC" and zincblende == True and surface == True:
                file_fdf.write("3  ")
            elif self.finalcell[i][3] == "H" and zincblende == False:
                file_fdf.write("2  ")
            elif self.finalcell[i][3] == "H" and zincblende == True and surface == False:
                file_fdf.write("3  ")
            elif self.finalcell[i][3] == "H" and zincblende == True and surface == True:
                file_fdf.write("4  ")
            # Below: Assignment of atoms relative to the other atoms and their coordinates
            if i == 0:
                file_fdf.write("0 0 0 ") # First atom therefore with respect nothing
                file_fdf.write(f"0.0 0.0 0.0 ") # first atom therefore with respect to itself
            elif i == 1:
                file_fdf.write("1 0 0 ") # with respect to the first atom
                polar = cartesian_to_polar(self.finalcell[i],self.finalcell[i-1])
                file_fdf.write(f"{polar[0]} {polar[1]} {polar[2]} ") # first atom therefore with respect to itself
            else:
                pass

            # Below: Constraining atoms
            if  self.finalcell[i][3] == "H":
                file_fdf.write("1 1 1")
            else:
                file_fdf.write("0 0 0")
            file_fdf.write("\n")
        file_fdf.write("%endblock Zmatrix\n")
        print("write_to_fdf_zmatrix: Successfully written to fdf file (z matrix)")
        file_fdf.close()

    def write_to_xmol(self, zincBlende, surface):
        # Printing the final into an out file that contains the coordinates
        file_xmol = open("coordinates.xyz", "w+")
        file_xmol.write(str(len(self.finalcell)) + "\n")
        file_xmol.write("The coordinates for the quantum dot atoms\n")
        if zincBlende:
            if surface:
                for i in self.finalcell:
                    if i[3] == "A":
                        file_xmol.write("A ")
                    elif i[3] == "C":
                        file_xmol.write("C ")
                    elif i[3] == "SA" or i[3] == "SC":
                        file_xmol.write("S ")
                    elif i[3] == "H":
                        file_xmol.write("H ")
                    file_xmol.write(str(i[0]) + " ")
                    file_xmol.write(str(i[1]) + " ")
                    file_xmol.write(str(i[2]) + " ")
                    file_xmol.write("\n")
            if not surface:
                for i in self.finalcell:
                    if i[3] == "A":
                        file_xmol.write("A ")
                    elif i[3] == "C":
                        file_xmol.write("C ")
                    elif i[3] == "SA":
                        file_xmol.write("A ")
                    elif i[3] == "SC":
                        file_xmol.write("C ")
                    elif i[3] == "H":
                        file_xmol.write("H ")
                    file_xmol.write(str(i[0]) + " ")
                    file_xmol.write(str(i[1]) + " ")
                    file_xmol.write(str(i[2]) + " ")
                    file_xmol.write("\n")

        elif not zincBlende:
            if surface:
                for i in self.finalcell:
                    if i[3] == "A":
                        file_xmol.write("A ")
                    elif i[3] == "C":
                        file_xmol.write("A ")
                    elif i[3] == "SA" or i[3] == "SC":
                        file_xmol.write("S ")
                    elif i[3] == "H":
                        file_xmol.write("H ")
                    file_xmol.write(str(i[0]) + " ")
                    file_xmol.write(str(i[1]) + " ")
                    file_xmol.write(str(i[2]) + " ")
                    file_xmol.write("\n")
            if not surface:
                for i in self.finalcell:
                    if i[3] == "A":
                        file_xmol.write("A ")
                    elif i[3] == "C":
                        file_xmol.write("A ")
                    elif i[3] == "SA":
                        file_xmol.write("A ")
                    elif i[3] == "SC":
                        file_xmol.write("A ")
                    elif i[3] == "H":
                        file_xmol.write("H ")
                    file_xmol.write(str(i[0]) + " ")
                    file_xmol.write(str(i[1]) + " ")
                    file_xmol.write(str(i[2]) + " ")
                    file_xmol.write("\n")
        print("write_to_xmol: Successfully written to xmol file")
        file_xmol.close()

    def write_to_nn(self):
        # Printing the final into an out file that contains the coordinates
        file_nn = open("coordinates.nn", "w+")
        try:
            for atom1 in self.finalcell:
                file_nn.write(
                    str(atom1[0]) + "  " + str(atom1[1]) + "  " + str(atom1[2]) + "  " + str(atom1[4]) + "\n")
                for atom2 in self.finalcell:
                    vec = abs(atom1[0]-atom2[0]) ** 2 + abs(atom1[1]-atom2[1]) ** 2 + abs(atom1[2]-atom2[2]) ** 2
                    if vec == 0.1875:
                        file_nn.write(
                            "nn:  " + str(atom2[0]) + "  " + str(atom2[1]) + "  " + str(atom2[2]) + "  " + str(
                                atom2[4]) + "\n")
                for atom2 in self.finalcell:
                    vec = abs(atom1[0]-atom2[0]) ** 2 + abs(atom1[1]-atom2[1]) ** 2 + abs(atom1[2]-atom2[2]) ** 2
                    if vec == 0.5:
                        file_nn.write(
                            "nnn:  " + str(atom2[0]) + "  " + str(atom2[1]) + "  " + str(atom2[2]) + "  " + str(
                                atom2[4]) + "\n")
            print("write_to_nn: Successfully written to nearest neighbour file")
        except:
            print("Warning: Could not write nearest neighbours. Please create atom numbers using the 'number_atoms' method")
        file_nn.close()

    def identify_site(self, *args):
        """This functions identifies a 'random' site that is of the specified arg types. List args by prefered order"""
        for arg in args:
            for x in self.finalcell:
                if x[3] == arg:
                    x[3] = "L"
                    print(f"idenify_site: Site found which is '{arg}'")
                    return x

    def attach_ligands(self, ligand):
        """Sites have to be identified for this function to work"""
        for x in self.finalcell:
            if x[3] == "L":
                pass
                # self.finalcell.append(x[3],[],[])



if __name__ == "__main__": 
    dot = DotDiamond(2, "Ang", 6.432, replicate=True)            # (Radius to build, lattice constant, replicate to all octants)
    dot.trim_to_dot(1.5, True)                                 # (radius to cutoff, evenize)
    # dot.oxygenate(1.7)                                           # Hydrogen bond length in Angstroms
    # dot.hydrogenate(1.7)                                         # Hydrogen bond length in Angstroms
    # print(dot.identify_site("H"))                                 # Coordinates of the site with 'H'
    # dot2_finalcell = invert_coordinates(dot.finalcell)           # Inverts the coordinates

# Do these things the last
    dot.number_atoms("creation")
    dot.write_to_log()
    dot.write_to_fdf(True, False)                               # (zincBlende, surface)
    # dot.write_to_fdf_zmatrix(False, False)                       # (zincBlende, surface)
    dot.write_to_xmol(True, False)                               # (zincBlende, surface)
    dot.write_to_nn()
    # print(f"The cartesian to polar conversion: {cartesian_to_polar([-1,-1,1],[0,0,0])}")
    # print(f"The angle between those three atoms: {triatom_angle([0,0,0],[1,1,0],[0,0,1])}")