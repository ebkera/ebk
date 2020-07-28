from datetime import datetime
import numpy as np
from numpy import linalg as LA

class TB:
    def __init__(self, *args, **kwargs):
        """
        atoms = variable stores all the atoms
        and Example for H_0
              --                    --
              | E_0     H_12    0    |
        H_0 = | H_12^T  E_0     H_12 |
              | 0       H_12^T  E_0  |
              --                    --
        """
        self.run_name = kwargs.get("run_name", "RUN")
        self.atoms = []
        self.so = kwargs.get("so", True)
        self.a = kwargs.get("a", 6.453)  # Lattice parameter in Angstroms
        fractional_bondlength = 0.25
        self.cations_as_neighbours = np.array([[fractional_bondlength,fractional_bondlength,fractional_bondlength],[fractional_bondlength,-fractional_bondlength,-fractional_bondlength],[-fractional_bondlength,fractional_bondlength,-fractional_bondlength],[-fractional_bondlength,-fractional_bondlength,fractional_bondlength]])*self.a
        self.nn_threshold = self.a*np.sqrt(3)/4  # the threshold distance for Nearest neighbours
        self.print_hamiltonian = False

        log_open = open(f"{self.run_name}.out", "w")
        log_open.write(f"#######################\n#### TIGHT BINDING ####\n##### Calculation  ####\n#######################\n{datetime.now()}\n\n")
        self.start_time = datetime.now()
        log_open.close()

    def log(self, text):
        log_open = open(f"{self.run_name}.out", "a")
        log_open.write(text)
        log_open.write("\n")
        log_open.close()

    def read_coordinates(self, file_name):
        from ase.io import read, write
        # file = ase.io.read(f"{path}", format = "espresso-out")
        # file = ase.io.read(f"{path}.out", format = "espresso-out")
        if "xyz" in file_name:
            with open(file_name, "r") as coordinate_file:
                for line in coordinate_file:
                    if len(line.split()) == 4:
                        temp = line.split()
                        self.atoms.append([temp[0], float(temp[1]), float(temp[2]), float(temp[3])])
        self.log(f"Coordinates have been read from file: {file_name}")
        self.log(f"------------ Reprinting coordinates and atoms:")
        for x in self.atoms:
            self.log(f"{x}")
        self.log(f"------------ End of coordinates and atoms\n")
        self.N = list(range(len(self.atoms)))   # Here are the number of atoms
        self.S = [0,1]                     # Here are the spins
        self.O = list(range(10))           # Here are the number of orbitals
        self.H = [0]

    def check(self):
        return False

        # if "STRUCT_IN" in file_name:
            