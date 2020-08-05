"""This file contails some utilities to be used with SIESTA"""

def xyz2fdf(file_name, format):
    """
    This function takes a .xyz file and converts it into a fdf compliant format file
    |Inputs: file_name :(string) (without extension)
    |        format: (string) "Ang" if anstroms
    |output: file_name.fdf file
    """

    # with open(f"{file_name}.xyz", 'r') as file: # Use file to refer to the file object
    #     # data = file.read()

    file = open(f"{file_name}.xyz", 'r')
    data = [line for line in file]
    # print(data)
    atom_number = data[0].strip()
    comment = data[1].strip("\n")
    print(f"xyz2fdf:Number of atoms is: {atom_number}")
    print(f"xyz2fdf:Comment           : {comment}")
    file.close()

    data_to_write = []
    species = []
    for x in range(2,len(data)):
        y = data[x].split()
        if y[0] not in species:
            species.append(y[0])
        species_number = species.index(y[0])
        # print(species_number)
        data_to_write.append(f"{y[1]}   {y[2]}   {y[3]}   {species_number + 1}")

    if format == "Ang" or "ang" or "angstroms" or "Angstroms":
        format = "Ang"

    file = open(f"{file_name}.fdf", 'w')
    file.write("# Generated using the xyz2fdf utility by Eranjan in ebk.SIESTA\n")
    file.write(f"NumberOfAtoms    {atom_number}\n")
    file.write(f"AtomicCoordinatesFormat  {format}\n")
    file.write(f"%block AtomicCoordinatesAndAtomicSpecies\n")
    for line in data_to_write:
        file.write(f"{line}\n")
    file.write(f"%endblock AtomicCoordinatesAndAtomicSpecies")
    file.close()


def siesta_convergence_checker(file_name):
    import matplotlib.pyplot as plt
    file = open(f"{file_name}", 'r')
    data = [line for line in file]
    file.close()
    # print(data)
    iteration_number = []
    Eharris = []
    E_KS = []
    FreeEng = []
    dDmax = []
    Ef = []
    dHmax = []
    for line in data:
        if "scf" in line and "compute" not in line and "siesta" not in line and "Eharris" not in line and "Vacuum" not in line and "dfscf" not in line:
            # print(line)
            vals = line.split()
            iteration_number.append(int(vals[1]))
            Eharris.append(float(vals[2]))
            E_KS.append(float(vals[3]))
            FreeEng.append(float(vals[4]))
            dDmax.append(float(vals[5]))
            Ef.append(float(vals[6]))
            dHmax.append(float(vals[7]))

    # Eharris = [-x for x in Eharris]
    # E_KS = [-x for x in E_KS]
    # FreeEng = [-x for x in FreeEng]
    # # print(FreeEng)

    plt.plot(iteration_number, Eharris, label=("Harris"))
    plt.plot(iteration_number, E_KS, label=("Khon-Sham"))
    plt.plot(iteration_number, FreeEng, label=("Free"))
    plt.title("SCF Convergence")
    plt.xlabel("Iteration")
    plt.ylabel("Energy")
    plt.legend()
    plt.savefig(f"{file_name}_SCF_convergence.pdf")
    plt.show()

    plt.figure()
    plt.plot(iteration_number, dDmax, label=("dDmax"))
    plt.plot(iteration_number, Ef, label=("E$_f$"))
    plt.plot(iteration_number, dHmax, label=("dHmax"))
    plt.title("Convergence Parameters and Fermi level")
    plt.xlabel("Iteration")
    plt.ylabel("Energy")
    plt.legend()
    plt.savefig(f"{file_name}_Convergence_parameters.pdf")
    plt.show()

def struct2xyz(file_name):
    """
    This function converts a SIESTA struct out file to n xyz file.
    """
    from ase.io import read, write
    atoms = read(file_name)
    print(f"struct2xyz: STRUCT_OUT file imported: {atoms}")
    file_name = file_name.split(".")
    del file_name[-1]
    file_write_name = ".".join(file_name)
    write(f"{file_write_name}.xyz", atoms)
    print(f"struct2xyz: Written to .xyz file")
