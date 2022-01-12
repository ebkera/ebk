"""
This file reads the the file out_file_name.out and extracts/calculates data/values from it
"""

class SiestaReadOut():
    def __init__(self, out_file_name):
        self.out_file_name = out_file_name
        self.out_file = []
        self.EIG_file = []
        self.Species = [["Label", "Atomic_Number_int", "Number_of_atoms_int"]]  # Logs all labels and atomic numbers. Since we have a example value at index=0 Index will serve as species number
        
        # Eig file stuff
        self.Ef_EIG = 0   # Fermi level as read in from eig
        self.eigen_values_by_k_point = [[]]  # Evey k point will have a set of eigen values. We have already filled one here for the 0th k point so from now index_of_list_element=kpoint

        # Here are defaults
        trigger_read_in_atoms = False

        f = open(f"{self.out_file_name}.out", "r")
        for line in f:
            self.out_file.append(line)
        f.close()
        f = open(f"{self.out_file_name}.EIG", "r")
        for line in f:
            self.EIG_file.append(line)
        f.close()

        # We will be reading all values here...
        for line in self.out_file:
            #Getting the number of types of atoms and masses from line example "Species number:   1 Atomic number:   50 Label: Sn"
            if "SystemLabel" in line:
                x = line.split()
                self.SystemLabel = x[1]
            if "Species number:" in line and "Atomic number:" in line:
                temp_Species = []
                x = line.strip("Species number:")
                x = x.split(":")
                x[0] = x[0].strip("Atomic number")
                x[1] = int(x[1].strip("Label"))
                x[2] = x[2].strip()
                self.Species.append([x[2], x[1], 0])
                # print(x)

            # Here we will read in vacuum levels
            if "dhscf: Vacuum level (max, mean)" in line:
                x = line.strip("dhscf: Vacuum level (max, mean)=")
                x = x.split()
                self.Vac_max = float(x[0])
                self.Vac_mean = float(x[1])
                self.Vac_units = x[2]

            # Reading fermi level here
            if "siesta:         Fermi =" in line:
                x = line.strip("siesta:         Fermi =")
                x = x.split()
                self.Ef = float(x[0])

            # Here we are reading in the Total energies
            if "siesta:         Total =" in line:
                x = line.strip("siesta:         Total =")
                x = x.split()
                self.E_total = float(x[0])      

            # Here we are reading in the number of atoms and orbitals and projectors
            if "initatomlists: Number of atoms, orbitals, and projectors:" in line:
                x = line.strip("initatomlists: Number of atoms, orbitals, and projectors:")
                x = x.split()
                self.N_atoms = int(x[0])
                self.N_orbitals = int(x[1])
                self.N_projectors = int(x[2])

            # Here we read in atoms should be below the actual reading in script
            if trigger_read_in_atoms:
                x = line.split()
                # print(x)
                if len(x) != 6: 
                    trigger_read_in_atoms = False
                    continue
                self.Species[int(x[4])][2] += 1


            # Here we set trigger to read in atoms (initial coordinates) should be below the actual reading in script
            if "siesta: Atomic coordinates (Bohr) and species" in line: trigger_read_in_atoms = True
            

# Here we are parsing in the .EIG file
        for i,line in enumerate(self.EIG_file):
            if i == 0:
                self.Ef_EIG = float(line)
                fermi_diff = abs(self.Ef - self.Ef_EIG)
                if fermi_diff >= 0.001:
                    print(f"Warning!: Fermi levels in the out file and EIG files don't match |{self.Ef_EIG} - {self.Ef}| = {fermi_diff}")
            if i == 1:
                x = line.split()
                self.NumberOfBands = int(x[0])
                self.NumberOfSpins = int(x[1])
                self.NumberOfkPoints = int(x[2])
            if i > 1:
                x = line.split()
                if len(x) == 11:
                    # This is a fresh K point line
                    new_eig_vals = []
                    for j,v in enumerate(x):
                        if j > 0:
                            new_eig_vals.append(float(v))                           
                else:
                    for j,v in enumerate(x):
                        new_eig_vals.append(float(v))
                self.eigen_values_by_k_point.append(new_eig_vals) 

        # While we are at it lets calculate the HOMO and LUMO levels. Not intelligent enough to identify conductors so keep in mind
        for i,x in enumerate(self.eigen_values_by_k_point[1]):
            if x >= self.Ef:
                self.band_index_LUMO = i+1       # Adding 1 to the index so that we the index starts from 1
                self.band_index_HOMO = (i-1)+1   # Adding 1 to the index so that we the index starts from 1
                break

                                        

# From here on we are doing methods
        # Here we are getting the Work functions
        try: 
            self.WF = self.Vac_mean - self.Ef
        except: print("SiestaReadOut: No vacuum values. Possible reasons: bulk system")
        # print(self.N_atoms, self.N_orbitals, self.N_projectors)
        # print(self.Species)

    def get_vacuum(self):
        """returns [self.Vac_max, self.Vac_mean, self.Vac_units]""" 
        return [self.Vac_max, self.Vac_mean, self.Vac_units]

    def get_total_energy(self):
        """returns self.E_total"""
        return self.E_total

    def get_fermi(self):
        """returns self.Ef"""
        return self.Ef

    def get_number_of_atoms(self, label = ""):
        """
        Returns the number of atoms of 'label' type if 'label' is null returns total number of atoms
        """
        atom_type = ""
        atom_type = label
        if atom_type == "": return self.N_atoms
        for x in range(len(self.Species)):
            if self.Species[x][0] == atom_type: return self.Species[x][2]
        # result = filter(lambda x: x[0] == atom_type, self.Species)
        
    def get_number_of_orbitals(self):
        """return self.N_orbitals"""
        return self.N_orbitals

    def get_number_of_projectors(self):
        """return self.N_projectors"""
        return self.N_projectors

    def get_number_of_bands(self):
        """return self.NumerOfBands"""
        return self.NumberOfBands

    def get_number_of_spins(self):
        """return self.NumberOfSpins"""
        return self.NumberOfSpins

    def get_number_of_k_points(self):
        """return self.NumberOfkPoints"""
        return self.NumberOfkPoints

    def get_HOMO_band_index(self):
        """Calculates the index of the HOMO level"""
        return self.band_index_HOMO

    def get_LUMO_band_index(self):
        """Calculates the index of the LUMO level"""
        return self.band_index_LUMO


# From here below we have deprecated code kept for back-compatibility
    def read_vacuum(self):
        """
        |Deprecated
        |
        |This function sets the maximum and the mean values of the vacuum and the units. 
        |Inputs:
        |    None
        |Outputs:
        |    None
        """
        for line in self.file:
            if "dhscf: Vacuum level (max, mean)" in line:
                x = line.strip("dhscf: Vacuum level (max, mean)=")
                x = x.split()
                self.Vac_max = float(x[0])
                self.Vac_mean = float(x[1])
                self.Vac_units = x[2]
                break

    def read_fermi(self):
        """
        |Deprecated
        |
        |This function reads the fermi level. For now it is set only to read in eV
        |Inputs:
        |    None
        |Outputs:
        |    Float: The fermi energy in eV
        """
        self.Ef = 0
        for line in self.file:
            if "siesta:         Fermi =" in line:
                x = line.strip("siesta:         Fermi =")
                x = x.split()
                self.Ef = float(x[0])
                return self.Ef

    def calculate_work_function(self):
        """
        |Deprecated
        |        
        |This function calculates the work functions and saves in self.WF 
        |Inputs:
        |    None
        |Outputs:
        |    None
        """
        try:
            # If fermi and vacuum levels are calculated we can put just use them
            self.WF = self.Vac_mean - self.Ef
        except:
            # If fermi and vacuum levels are not calculated we need to calculate them first
            self.read_fermi()
            self.read_vacuum()
            self.WF = self.Vac_mean - self.Ef
        return self.WF

    def read_total_energy(self):
        """
        |Deprecated
        |
        |This function reads the total_energy. For now it is set only to read in eV
        |Inputs:
        |    None
        |Outputs:
        |    None
        """
        # self.Ef = 0
        for line in self.file:
            if "siesta:         Total =" in line:
                x = line.strip("siesta:         Total =")
                x = x.split()
                self.E_total = float(x[0])
                return self.E_total

def get_quadrupole_moments(file_name):
    """
    Calculates the quadrupole for the sytem given.
    Reads in the *.RHO.cube file generated by Denchar.


    The basic structure of the cube file is given here.
    http://paulbourke.net/dataformats/cube/

    Header
    The first two lines of the header are comments, they are generally ignored by parsing packages or used as two default labels.

    The third line has the number of atoms included in the file followed by the position of the origin of the volumetric data.

    The next three lines give the number of voxels along each axis (x, y, z) followed by the axis vector. Note this means the volume need not be aligned with the coordinate axis, indeed it also means it may be sheared although most volumetric packages won't support that. The length of each vector is the length of the side of the voxel thus allowing non cubic volumes. If the sign of the number of voxels in a dimension is positive then the units are Bohr, if negative then Angstroms.

    The last section in the header is one line for each atom consisting of 5 numbers, the first is the atom number, the second is the charge, and the last three are the x,y,z coordinates of the atom center.

    Volumetric data
    The volumetric data is straightforward, one floating point number for each volumetric element. The original Gaussian format arranged the values in the format shown below in the example, most parsing programs can read any white space separated format. Traditionally the grid is arranged with the x axis as the outer loop and the z axis as the inner loop, for example, written as

    for (ix=0;ix<NX;ix++) {
        for (iy=0;iy<NY;iy++) {
            for (iz=0;iz<NZ;iz++) {
                printf("%g ",data[ix][iy][iz]);
                if (iz % 6 == 5)
                printf("\n");
            }
            printf("\n");
        }
    }

    """
    import numpy as np

    read_file_name = file_name
    origin = []
    number_of_atoms = 0
    x_voxel_counter = 0
    y_voxel_counter = 0
    z_voxel_counter = 0
    x_number_of_voxels = 0
    y_number_of_voxels = 0
    z_number_of_voxels = 0
    x_voxel_length = 0
    y_voxel_length = 0
    z_voxel_length = 0
    atom_coordinates_tigger = False
    total_charge = 0

    with open(read_file_name, "r") as file:
        for i,line in enumerate(file):
            if i == 2:
                parsed = line.split()
                number_of_atoms = int(parsed[0])
                for x in range(1,len(parsed)):
                    origin.append(float(parsed[x]))
            if i == 3:
                parsed = line.split()
                x_number_of_voxels = int(parsed[0])
                x_voxel_length = np.sqrt(float(parsed[1])**2+float(parsed[1])**2+float(parsed[1])**2)
                # print(parsed)
            if i == 4:
                parsed = line.split()
                y_number_of_voxels = int(parsed[0])
                y_voxel_length = np.sqrt(float(parsed[1])**2+float(parsed[1])**2+float(parsed[1])**2)
                # print(parsed)
            if i == 5:
                parsed = line.split()
                z_number_of_voxels = int(parsed[0])
                z_voxel_length = np.sqrt(float(parsed[1])**2+float(parsed[1])**2+float(parsed[1])**2)
                # print(parsed)
                data = np.zeros((x_number_of_voxels, y_number_of_voxels, z_number_of_voxels))
                print(data.shape)

            
            if i == 6: atom_coordinates_tigger = True

            if atom_coordinates_tigger:
                if len(line.split()) != 5:
                    atom_coordinates_tigger = False

            if not atom_coordinates_tigger and i>6:
                parsed = line.split()
                for val in parsed:
                    data[x_voxel_counter,y_voxel_counter,z_voxel_counter] = float(val)
                    total_charge+= float(val)
                    z_voxel_counter+=1
                    # print(z_voxel_counter)
                    if z_voxel_counter == z_number_of_voxels: 
                        z_voxel_counter = 0
                        y_voxel_counter+=1
                        if y_voxel_counter == y_number_of_voxels:
                            y_voxel_counter=0
                            x_voxel_counter+=1


    print(data[12,13,14])
    print(total_charge)

    Qxx = 0
    Qyy = 0
    Qzz = 0

    for x in range(x_number_of_voxels):
        for y in range(y_number_of_voxels):
            for z in range(z_number_of_voxels):
                """Methana podi indeces proshnayak thiyeanwa. mokenda iterate karanna one i+1 da nattam i da kiyala"""
                r2 = (x*x_voxel_length+x_voxel_length/2)**2*(y*y_voxel_length+y_voxel_length/2)**2*(z*z_voxel_length+z_voxel_length/2)**2
                # Qxx+= data[x,y,z]*(3*x_voxel_length*x*x_voxel_length*x - r2)
                Qxx+= data[x,y,z]*(3*(x*x_voxel_length+x_voxel_length/2)*(x*x_voxel_length+x_voxel_length/2) - r2)
                Qyy+= data[x,y,z]*(3*(y*y_voxel_length+y_voxel_length/2)*(y*y_voxel_length+y_voxel_length/2) - r2)
                Qzz+= data[x,y,z]*(3*(z*z_voxel_length+z_voxel_length/2)*(z*z_voxel_length+z_voxel_length/2) - r2)


    print("Qxx", Qxx)
    print("Qyy", Qyy)
    print("Qzz", Qzz)


if __name__ == "__main__":
    vac = SiestaReadOut("Fe")
    # vac.read_vacuum()
    # vac.read_fermi()
    vac.read_total_energy()
    print(vac.E_total)

    vac.calculate_work_function()
    print(vac.Ef)
    print(vac.Vac_mean)
    print(vac.WF)
