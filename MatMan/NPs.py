import matplotlib.pyplot as plt
from ase.io import read, write
from ase import Atoms
from ebk.SIESTA import read_struct_file, struct2xyz
from ebk import get_machine_paths
import numpy as np
# from ebk.coordinateMaker import

def analyze_bonds(file_name, central_atom_index = 0, passivant = "H", dot_species = "Sn"):
    """
    Reads in an XYZ file.
    This function analyzes the following,
    The radis of the dot : this is crucial after relaxation especially.
    The relative bond lenghts.
    and Graphs the said information and also outputs a log file. with latex compatible stuff.
    structure of atoms list is [atom, radial_dist_from_central_atom, nn1, nn2, nn3, nn4, surface(bool)]
    """

    def write_to_log(fname = file_name):
        """Printing the final into an out file that contains the coordinates"""
        fname = fname.split(".")
        del fname[-1]
        fname = ".".join(fname)
        file_my = open(f"{fname}.log", "w+")        
        file_my.write("\n\n##--Log of run in Latex format--(Copy and paste this into .tex file)--\n\n")
        file_my.write("\\begin{table}[h]\n")
        file_my.write("\\centering\n")
        file_my.write("\\begin{tabular}{ll}\n")
        file_my.write("Coordinate Format             &: " + str(coordinate_format)+" \\\\\n")
        file_my.write("NP atom Count                 &: " + str(no_of_NP_atoms) + "\\\\\n")
        try:
            # file_my.write("Anion/Cation Ratio            &: " + str(anion_count/cation_count) + "\\\\\n")
            file_my.write("Anion Count                   &: " + str(anion_count) + "\\\\\n")
            file_my.write("Cation Count                  &: " + str(cation_count) + "\\\\\n")
        except: pass
        try:
            file_my.write("Surface sites                 &: " + str(no_of_surface_atoms) + "\\\\\n")
        except:
            file_my.write("% No Surface atoms detected!\n")
        try:
            file_my.write("Passivation atom count        &: " + str(no_of_passivant_atoms) + "\\\\\n")
            file_my.write("Passivation bond length       &: " + str(passivation_bondlength) + " (\\AA)\\\\\n")
        except:
            file_my.write("% The initial box was never pasivated! - No Passivation atoms!\n")
        try:
            file_my.write("Dot diameter                  &: " + str(diameter)+ " (\\AA) \\\\\n")
            file_my.write("Dot diameter                  &: " + str(diameter/10)+ " (nm) \\\\\n")
            # file_my.write("Dot diameter with passivation &: $\\approx$ " + str(diameter+passivation_bondlength)+ " (\\AA)\\\\\n")
        except:
            file_my.write("% The initial box was not trimmed! - No defined cutoff - No diameter or radius!\n")
        file_my.write("Total atoms                   &: " + str(len(atoms)) + "\n")
        file_my.write("\\end{tabular}\n")
        file_my.write("\\caption{QD with " + str(len(atoms)) + " total number of atoms} \\label{tab:QD" + str(len(atoms)) + "}\n")
        file_my.write("\\end{table}\n")
        file_my.close()
        print("write_to_log: Successfully written to the log file")

    # Reading in the file
    atoms_ase = read(file_name)
    coordinate_format = "Angstroms"
    atoms = []
    central_atom = atoms_ase[central_atom_index]
    bond_length_threshold = 1.1*atoms_ase.get_distance(0, 1)
    surface_radii = []
    no_of_NP_atoms = 0
    no_of_passivant_atoms = 0
    for i, atom in enumerate(atoms_ase):
        if atom.symbol == dot_species:
            no_of_NP_atoms+=1
        elif atom.symbol == passivant:
            no_of_passivant_atoms += 1
        list_to_append = []
        surface = False
        list_to_append.append(atom)
        list_to_append.append(atoms_ase.get_distance(0, i))
        # Looking for nearest neighbours
        for j, atom2 in enumerate(atoms_ase):
            if i != j and atoms_ase.get_distance(i, j) < bond_length_threshold:
                list_to_append.append(atom2)
                if atom2.symbol == "H": surface = True
        list_to_append.append(surface)
        atoms.append(list_to_append)

    # Getting Dot radius
    for atom in atoms:
        # print(atom[1])s
        if atom[0].symbol == "Sn":
            if atom[6] == True: surface_radii.append(atom[1])

    no_of_surface_atoms = len(surface_radii)
    radius = sum(surface_radii)/no_of_surface_atoms
    diameter = 2*radius
    print(f"Diameter: {diameter} (Angstroms)")
    write_to_log()

def create_shell_NP(file_name, central_atom_index = 0, passivant = "H", dot_species = "Sn", bond_length = False, N_shells = 1):
    """
    eads in an XYZ file.
    This function does the following,
    Creates a given number of shells for a given nano-particle
    """
    fname = file_name.split(".")
    del fname[-1]
    fname = ".".join(fname)
    
    # Reading in the file
    atoms_ase = read(file_name)
    coordinate_format = "Angstroms"
    atoms = []
    central_atom = atoms_ase[central_atom_index]
    if not bond_length: bond_length = atoms_ase.get_distance(0, 1)
    bond_length_threshold = 1.2*atoms_ase.get_distance(0, 1)
    surface_radii = []
    no_of_NP_atoms = 0
    no_of_passivant_atoms = 0
    surface_atom_indeces = []
    passivant_atom_indeces = []
    for i, atom in enumerate(atoms_ase):
        if atom.symbol == dot_species:
            no_of_NP_atoms+=1
        elif atom.symbol == passivant:
            no_of_passivant_atoms += 1
            passivant_atom_indeces.append(i)
            continue
        # list_to_append = []
        # surface = False
        # list_to_append.append(atom)
        # list_to_append.append(atoms_ase.get_distance(0, i))
        # Looking for nearest neighbours
        for j, atom2 in enumerate(atoms_ase):
            if i != j and atoms_ase.get_distance(i, j) < bond_length_threshold:
                # list_to_append.append(atom2)
                if atom2.symbol == "H":
                    if i not in surface_atom_indeces:   # to avoind appendign teh same surface atoms that is connected to multiple passivation atoms
                        surface_atom_indeces.append(i)
                    atom.symbol = "S"
        # list_to_append.append(surface)
        # atoms.append(list_to_append)

    # Printing out the NP with surface atoms highlighted
    write(f"{fname}_SA.xyz", atoms_ase)
    #removing the H atoms
    passivant_atom_indeces.sort(reverse = True) # sorting in reverse order in anticipation of deletion
    for i in passivant_atom_indeces:
        del(atoms_ase[i])
    write(f"{fname}_!PA.xyz", atoms_ase)

    # Now starts the layer growing stuff

    default_positions = [[bond_length,bond_length,bond_length],[-bond_length,-bond_length,bond_length],[-bond_length,bond_length,-bond_length],[bond_length,-bond_length,-bond_length]]
    added_atoms_list = []

    for layer in range(1,N_shells+1):
        print(f"Working on layer: {layer}")

        for i in surface_atom_indeces:
            # print(atoms_ase[i].position)
            radial_dist = atoms_ase.get_distance(0, i)
            for pos in default_positions:
                new_position = [atoms_ase[i].position[0] + pos[0], atoms_ase[i].position[1] + pos[1], atoms_ase[i].position[2] + pos[2]]
                new_dist = np.sqrt(new_position[0]**2+new_position[1]**2+new_position[2]**2)
                if new_dist < radial_dist: continue
                atoms_ase.append("C") # for newly added
                atoms_ase[-1].position = [atoms_ase[i].position[0] + pos[0], atoms_ase[i].position[1] + pos[1], atoms_ase[i].position[2] + pos[2]]
                added_atoms_list.append(len(atoms_ase)-1)
        print(added_atoms_list)

        # cleaning up the unwanted atoms
        atoms_to_delete = []
        for i in added_atoms_list:
            for j, atom2 in enumerate(atoms_ase):
                if i != j and atoms_ase.get_distance(i, j) < bond_length_threshold:
                    atoms_to_delete.append(j)

    write(f"{fname}_final.xyz", atoms_ase)


def create_relaxed(file_name, replace_with):
    """
    Reads in an XYZ file.
    This function does the following,
    Creates a given number of shells for a given nano-particle
    replacewith should also be a file name.
    """
    fname = file_name.split(".")
    del fname[-1]
    fname = ".".join(fname)
    
    # Reading in the file
    atoms_ase = read(file_name)
    replace_with = read(replace_with)

    bond_length_threshold = atoms_ase.get_distance(0, 1)

    list_to_delete = []
    for i, atom in enumerate(replace_with):
        if atom.symbol == "H":
            list_to_delete.append(i)
    list_to_delete.sort(reverse = True)
    for i in list_to_delete:
        del(replace_with[i])

    list_to_delete = []
    for i, atom in enumerate(replace_with):
        for j, atom2 in enumerate(atoms_ase):
            if atom2.symbol == "H": continue
            distx = abs(atom.position[0]-atom2.position[0])
            disty = abs(atom.position[1]-atom2.position[1])
            distz = abs(atom.position[2]-atom2.position[2])
            dist = np.sqrt(distx**2+disty**2+distz**2)
            if dist < 0.5*bond_length_threshold: 
                print(f"ase: {atom2} is close to replace atom {atom}")
                list_to_delete.append(j)

    list_to_delete.sort(reverse = True)
    for i in list_to_delete:
        del(atoms_ase[i])

    atoms_ase.extend(replace_with)

    write(f"{fname}_SemiRelaxed.xyz", atoms_ase)




class Insert_NP():
    """
    This class is made to be strictly for ligands so that we dont have to go about attaching ligands manually. 
    This class will insert a ligand at a given position and at a given orientation.
    """
    def __init__(self, *args, **kwargs):
        self.species = kwargs.get("species", "Sn")
        self.size = kwargs.get("size", 47)
        self.name = kwargs.get("name", f"{self.species}{self.size}")
        self.relaxed = kwargs.get("relaxed", True)
        self.XC = kwargs.get("XC", "LDA")
        self.SO_strength = kwargs.get("SO_strength", 1)

        # Loading the Nano Particle
        path = get_machine_paths()["xyz"]

        if self.relaxed == True:
            self.atoms = read_struct_file(f"{path}/NPs/Rivero/LDA_Relaxed_nonriveroPAObasis_SO_1/{self.species}{self.size}_{self.XC}_SOstrength_{self.SO_strength}.STRUCT_OUT")

    def orient(self, direction):
        """
        Plan to make this into a function that takes in  direction and orients the ligands in that direction.
        This might have to be a ligands specific fucntion though.
        """
        if direction == [0,0,1]:
            self.atoms.rotate(90, 'y')
        if direction == [1,1,0]:
            self.atoms.rotate(45, 'z')

    def edit(self):
        self.atoms.edit()

    def get_chemical_symbols(self):
        return self.atoms.get_chemical_symbols()

    def get_positions(self):
        return self.atoms.get_positions()

    def list_all_atoms(self):
        syms = self.get_chemical_symbols()
        coor = self.get_positions()
        return [f"{syms[n]}: {coor[n]}" for n in range(len(syms))]

    def save(self, format):
        write(f"{self.name}.{format}", self.atoms)

    def attach(self, attach_to, attach_at, attach_through):
        """
        This is the main function that attacheds the ligand to a atoms type object
        attach_to: The atoms type object that the ligand will attach to
        attach_at: Will attach at this atom in the atoms object(attach_to) and if there is an atom there already it will remove it and attach the ligand
        attach_through: The ligand will remove this atom and the resulting will attach throught the resulting dangling bond. to the attach_at site in the ligand.
        """
        pass

class NP(Insert_NP):
    def __init__(self, *args, **kwargs):
        super().__init__(self, *args, **kwargs)

# class BDT12(Insert_ligand):
#     def __init__(self, *args, **kwargs):
#         self.name = kwargs.get("name", "BDT12")
#         self.atoms = read(BDT12_path, index=None, format="xyz")
#         super().__init__()

# class EDT12(Insert_ligand):
#     def __init__(self, *args, **kwargs):
#         self.name = kwargs.get("name", "EDT12")
#         self.atoms = read(EDT12_path, index=None, format="xyz")
#         super().__init__()