"""This file contails some utilities to be used with SIESTA"""

def xyz2fdf(file_name, format, lattice=False, lattice_constant=0):
    """
    This function takes a .xyz file and converts it into a fdf compliant format file
    |Inputs: file_name :(string) (without extension)
    |        format: (string) "Ang" if anstroms
    |        lattice: (bool) If true will try to write the lattice vectors from data in the xyz file.
    |        lattice constant: (float) Needed if Lattice is True.
    |output: file_name.fdf file
    """

    # with open(f"{file_name}.xyz", 'r') as file: # Use file to refer to the file object
    #     # data = file.read()

    if lattice and lattice_constant == 0:
        raise Exception("Sorry Lattice Constant not set but lattice vectors are expected to be written. This will make SIESTA treat this as a molecule. Please set lattice_constant or set lattice=False. ")

    file = open(f"{file_name}.xyz", 'r')
    data = [line for line in file]
    # print(data)
    atom_number = data[0].strip()
    comment = data[1].strip("\n")
    print(f"xyz2fdf:Number of atoms is: {atom_number}")
    print(f"xyz2fdf:Comment           : {comment}")
    file.close()

    if "Lattice" in comment and lattice:
        # import re
        # txt = re.sub("[0-9]+", "x", comment)
        # txt = re.findall(r'[0-9]+', comment) 
        txt = comment.strip("Lattice")
        txt = txt.strip('="')
        txt = txt.split(" ")
        a = []
        b = []
        c = []
        for x in range(0,9):
            y = float(txt[x].strip('"'))/lattice_constant
            if x < 3: a.append(y)
            elif x<6: b.append(y)
            elif x<9:c.append(y)
    data_to_write = []
    species = []
    for x in range(2,len(data)):
        y = data[x].split()
        if y[0] not in species:
            species.append(y[0])
        species_number = species.index(y[0])
        # print(species_number)
        data_to_write.append(f"{float(y[1]):3.12f}\t{float(y[2]):3.12f}\t{float(y[3]):3.12f}\t{species_number + 1:>}")

    if format == "Ang" or "ang" or "angstroms" or "Angstroms":
        format = "Ang"

    file = open(f"{file_name}.fdf", 'w')
    file.write("# Generated using the xyz2fdf utility by Eranjan in ebk.SIESTA\n")
    file.write(f"NumberOfAtoms    {atom_number}\n")
    file.write(f"AtomicCoordinatesFormat  {format}\n")

    if lattice:
        file.write(f"LatticeConstant          {lattice_constant} Ang\n\n")
        if 'Lattice' in comment:
            file.write(f"%block LatticeVectors\n")
            file.write(f"{a[0]:.12f}  {a[1]:.12f}  {a[2]:>.12f}\n")
            file.write(f"{b[0]:.12f}  {b[1]:.12f}  {b[2]:>.12f}\n")
            file.write(f"{c[0]:.12f}  {c[1]:.12f}  {c[2]:>.12f}\n")
            file.write(f"%endblock LatticeVectors\n\n")
        else:
            file.write(f"#  xyz2fdf: WARNING: Tried to read lattice from xyz comments but failed!!!\n\n")
            
    file.write(f"%block AtomicCoordinatesAndAtomicSpecies\n")
    for line in data_to_write:
        file.write(f"{line}\n")
    file.write(f"%endblock AtomicCoordinatesAndAtomicSpecies\n\n")
    file.close()

def siesta_convergence_checker(file_name, title_addon="",show_linear_Kicks=False, show_struct_opt_moves=False, show_parameters=False, show_Harris=False):
    """
    This is legacy code now
    Use convergence_checker method in the siesta SIESTAOutFileReader.py    
    """
    import matplotlib.pyplot as plt
    file = open(f"{file_name}", 'r')
    data = [line for line in file]
    file.close()

    text_to_write = "# Eranjan\n"
 
    # print(data)
    iteration_number = []
    inter_num_count = 1
    scf_num = []
    Eharris = []
    E_KS = []
    FreeEng = []
    dDmax = []
    Ef = []
    dHmax = []
    Liner_kick_at = []
    Struct_opt_moves = []
    for line in data:
        if "scf" in line and "compute" not in line and "siesta" not in line and "Eharris" not in line and "Vacuum" not in line and "dfscf" not in line and "spin moment" not in line:
            # print(line)
            try:
                vals = line.split()
                if len(vals) != 8 : continue
                dDmax.append(float(vals[5]))  # This is here uptop because somtime it breaks for MD steps here and will go into the except before iterating iteration_number
                iteration_number.append(inter_num_count)
                scf_num.append(int(vals[1]))
                inter_num_count+=1
                Eharris.append(float(vals[2]))
                E_KS.append(float(vals[3]))
                FreeEng.append(float(vals[4]))
                Ef.append(float(vals[6]))
                dHmax.append(float(vals[7]))
            except: continue
            text_to_write+=line
            text_to_write+=f"len_vals|len_iter|len_dDmax|len_Ef:{len(vals)|len(iteration_number)}|{len(dDmax)}|{len(Ef)}\n"
        if "Linear-Kick" in line and "switching mixer" in line:
            text_to_write+=line
            text_to_write+=f"linear kick here ->\n"
            Liner_kick_at.append(inter_num_count)
        if "Begin" in line and "opt. move =" in line:
            text_to_write+=line
            text_to_write+=f"Opt Move here ->\n"
            if inter_num_count == 1: Struct_opt_moves.append(1)
            else: Struct_opt_moves.append(inter_num_count+1)

    file = open(f"{file_name}_scf_convergence.era", 'w')
    file.write(text_to_write)
    file.close()

    if show_Harris: plt.plot(iteration_number, Eharris, 'c', label=("Harris"))
    plt.rcParams["figure.figsize"] = (120,12)
    plt.plot(iteration_number, FreeEng, 'g', label=("Free Energy"))
    plt.plot(iteration_number, E_KS, "b",label=("Kohn-Sham"))
    plt.title(f"SCF Convergence: {title_addon}")
    plt.xlabel("Iteration")
    plt.ylabel("Energy (eV)")
    if show_linear_Kicks:
        if len(Liner_kick_at) != 0:
            plt.axvline(Liner_kick_at[0], color='r', linestyle='-', linewidth=.5, label="Linear Kicks")
            for x in range(1,len(Liner_kick_at)):
                plt.axvline(Liner_kick_at[x], color='r', linestyle='-', linewidth=.5, )
    if show_struct_opt_moves:
        if len(Struct_opt_moves) != 0:
            plt.axvline(Struct_opt_moves[0], color='m', linestyle='-', linewidth=.5, label="Opt. move")
            for x in range(1,len(Struct_opt_moves)):
                plt.axvline(Struct_opt_moves[x], color='m', linewidth=.5, linestyle='-')                
    plt.legend()
    plt.savefig(f"{file_name}_SCF_convergence.pdf")
    plt.show()

    if show_parameters:
        plt.figure()
        if show_linear_Kicks:
            if len(Liner_kick_at) != 0:
                plt.axvline(Liner_kick_at[0], color='r', linestyle='-', linewidth=.5, label="Linear Kicks")
                for x in range(1,len(Liner_kick_at)):
                    plt.axvline(Liner_kick_at[x], color='r', linestyle='-', linewidth=.5, )
        if show_struct_opt_moves:
            if len(Struct_opt_moves) != 0:
                plt.axvline(Struct_opt_moves[0], color='m', linestyle='-', linewidth=.5, label="Opt. move")
                for x in range(1,len(Struct_opt_moves)):
                    plt.axvline(Struct_opt_moves[x], color='m', linewidth=.5, linestyle='-')                
        plt.plot(iteration_number, dDmax, label=("dDmax"))
        plt.plot(iteration_number, Ef, label=("E$_f$"))
        plt.plot(iteration_number, dHmax, label=("dHmax"))
        plt.title(f"Convergence Parameters and Fermi level: {title_addon}")
        plt.xlabel("Iteration")
        plt.ylabel("Energy (eV)")
        plt.legend()
        plt.savefig(f"{file_name}_Convergence_parameters.pdf")
        plt.show()

def struct2xyz(file_name):
    """
    This function converts a SIESTA struct out file to an xyz file.
    Inputs: filename (string): The file name to save the file. Should the given without the extension.
    """
    from ase.io import read, write
    atoms = read(file_name, format="struct_out")
    print(f"struct2xyz: STRUCT_OUT file imported: {atoms}")
    file_name = file_name.split(".")
    del file_name[-1]
    file_write_name = ".".join(file_name)
    write(f"{file_write_name}.xyz", atoms)
    print(f"struct2xyz: Written to .xyz file")

def struct2vasp(file_name):
    """
    This function converts a SIESTA struct out file to a vasp file.
    Inputs: filename (string): The file name to save the file. Should the given without the extension.
    """
    from ase.io import read, write
    atoms = read(file_name, format="struct_out")
    print(f"struct2vasp: STRUCT_OUT file imported: {atoms}")
    file_name = file_name.split(".")
    del file_name[-1]
    file_write_name = ".".join(file_name)
    write(f"{file_write_name}.vasp", atoms)
    print(f"struct2vasp: Written to .vasp file")

def struct2cif(file_name):
    """
    This function converts a SIESTA struct out file to a vasp file.
    Inputs: filename (string): The file name to save the file. Should the given without the extension.
    """
    from ase.io import read, write
    atoms = read(file_name, format="struct_out")
    print(f"struct2vasp: STRUCT_OUT file imported: {atoms}")
    file_name = file_name.split(".")
    del file_name[-1]
    file_write_name = ".".join(file_name)
    write(f"{file_write_name}.cif", atoms)
    print(f"struct2vasp: Written to .cif file")

def read_struct_file(file_name):
    """
    This function reads a struct file and returns the atoms object.
    """
    from ase.io import read
    atoms = read(file_name, format="struct_out")
    print(f"read_struct_file: STRUCT_OUT file imported: {atoms}")
    return atoms

def get_geometrical_steps(file_name:str, write_png:bool = False, animate:bool=True) -> None:
    """
    This function breaks up *.ANI (.ANI files are in xyz format) files into multiple .xyz files with subscripts.
    file_name: (str) should be the file name of the .ANI file (or path) should include the .ANI
    write_png: (bool) should be True if every step is to be written to a png file.
    animate  : (bool) Setting this to true makes .mp4 and .gif animations
    """
    from ase.io import read, write
    from ase.io.animation import write_gif, write_mp4

    filename_parts = file_name.split(".")
    # Check to see if it is a .ANI file, not completed
    del filename_parts[-1]
    filename_pre = ".".join(filename_parts)
    from ase.io import read
    file_number = 0
    atom_images_z = []
    atom_images_x = []
    atom_images_y = []
    with open(file_name, "r+") as file:
        for line in file:
            len_of_line = len(line.split())
            if len_of_line == 1:
                try:
                    file_to_write.close()
                    if write_png:
                        atoms = read(f"{filename_pre}.{file_number}.xyz")
                        write(f"{filename_pre}.z{file_number}.png", atoms)
                        atoms.rotate([1,0,0], [0,0,1])
                        write(f"{filename_pre}.x{file_number}.png", atoms)
                        atoms.rotate([0,1,0], [0,0,1])
                        write(f"{filename_pre}.y{file_number}.png", atoms)
                    if animate:
                        atoms = read(f"{filename_pre}.{file_number}.xyz")
                        atom_images_z.append(atoms)
                        atomsx = atoms.copy()
                        atomsx.rotate([1,0,0], [0,0,1])
                        atom_images_x.append(atomsx)
                        atomsy = atomsx.copy()
                        atomsy.rotate([0,1,0], [0,0,1])
                        atom_images_y.append(atomsy)
                except:
                    pass
                file_number+=1
                file_to_write = open(f"{filename_pre}.{file_number}.xyz", "w+")
            file_to_write.write(line)
        file_to_write.close()
        print(f"Number of Frames: {len(atom_images_y)}")
        if write_png:
            atoms = read(f"{filename_pre}.{file_number}.xyz")
            write(f"{filename_pre}.z{file_number}.png", atoms)
            atoms.rotate([1,0,0], [0,0,1])
            write(f"{filename_pre}.x{file_number}.png", atoms)
            atoms.rotate([0,1,0], [0,0,1])
            write(f"{filename_pre}.y{file_number}.png", atoms)
        if animate:
            print(f"Making mp4 animations...")
            write_mp4(f"{filename_pre}.z.mp4", atom_images_z)
            write_mp4(f"{filename_pre}.x.mp4", atom_images_x)
            write_mp4(f"{filename_pre}.y.mp4", atom_images_y)
            print(f"Making gif animations...")
            write_gif(f"{filename_pre}.z.gif", atom_images_z)
            write_gif(f"{filename_pre}.x.gif", atom_images_x)
            write_gif(f"{filename_pre}.y.gif", atom_images_y)


def get_cwd():
    import os
    return os.getcwd()

def get_run_name():
    cwd = get_cwd()
    return cwd.split("/")[-1].split("^")[0]