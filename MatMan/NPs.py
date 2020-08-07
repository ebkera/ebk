import matplotlib.pyplot as plt
from ase.io import read, write
from ase import Atoms
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




    # for atom in atoms:
    #     print(len(atom))



    # 
