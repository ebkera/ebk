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
    KS_energy = []
    for line in data:
        if "scf" in line and "compute" not in line and "siesta" not in line and "Eharris" not in line:
            # print(line)
            vals = line.split()
            iteration_number.append(int(vals[1]))
            KS_energy.append(float(vals[3]))

    plt.plot(iteration_number, KS_energy)
    plt.title("SCF Convergence")
    plt.xlabel("Iteration")
    plt.ylabel("Energy")
    plt.show()
    
