"""
This file reads the the file out_file_name.out and extracts/calculates data/values from it
"""
from numpy import char
import scipy
from ebk import progress_bar
class SiestaReadOut():
    def __init__(self, out_file_name):
        self.out_file_name = out_file_name
        self.out_file = []
        self.EIG_file = []
        self.bands_file = []
        self.RHO_file = []
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
        try:
            f = open(f"{self.out_file_name}.bands", "r")
            for line in f:
                self.bands_file.append(line)
            f.close()
        except: pass

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

            # # Here we read in atoms should be below the actual reading in script
            # if trigger_read_in_atoms:
            #     x = line.split()
            #     print(x)
            #     if len(x) != 6: 
            #         trigger_read_in_atoms = False
            #         continue
            #     self.Species[int(x[4])][2] += 1
            # # Here we set trigger to read in atoms (initial coordinates) should be below the actual reading in script
            # if "siesta: Atomic coordinates (Bohr) and species" in line: trigger_read_in_atoms = True

            # Here we read in the number of electrons
            if "Total number of electrons" in line:
                self.number_of_electrons = float(line.split()[4])

            # Here we read in the total charge
            if "Total ionic charge" in line:
                self.total_ionic_charge = float(line.split()[3])
            

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

        # Parsing the bands files
        self.bands = []
        bands_temp = []
        self.k_dist = []
        self.hsp = []
        self.hss = []
        k_point = 0
        self.highest_valance = [0, -500]    # Will contain [kdist,E]
        self.lowest_conduction = [0, 500]   # Will contain [kdist,E]        
        for i,line in enumerate(self.bands_file):
            parsed = line.split()
            if i == 0:
                self.Ef_bands = float(parsed[0])
            if i == 1:
                self.kmin = float(parsed[0])
                self.kmax = float(parsed[1])
            if i == 2:
                self.Emin = float(parsed[0])
                self.Emax = float(parsed[1])
            if i == 3:
                self.NumberOfBands_bands = float(parsed[0])
                self.NumberOfSpins_bands = float(parsed[1])
                self.NumberOfkPoints_bands = float(parsed[2])  # this is different to self.NumberOfKPoints so we
            if i>3:
                if len(parsed) != 11 and len(bands_temp) == 0 and len(parsed) == 1:
                    self.NumberOfkLines = float(parsed[0])
                    continue
                if len(parsed) != 11 and len(bands_temp) == 0 and len(parsed) == 2:
                    self.hsp.append(float(parsed[0]))
                    parsed[1] = parsed[1].strip("'")
                    if parsed[1] == "GAMMA": parsed[1] = "$\Gamma$"
                    self.hss.append((parsed[1]))
                    continue
                if len(parsed) == 11:
                    self.k_dist.append(float(parsed.pop(0)))
                for val in parsed:
                    val = float(val) - self.Ef_bands
                    bands_temp.append(val)
                    if val < self.lowest_conduction[1] and len(bands_temp) == self.get_LUMO_band_index():
                            self.lowest_conduction[1] = val
                            self.lowest_conduction[0] = self.k_dist[-1]
                    if val > self.highest_valance[1] and len(bands_temp) == self.get_HOMO_band_index():
                            self.highest_valance[1] = val
                            self.highest_valance[0] = self.k_dist[-1]
                if len(bands_temp) == self.NumberOfBands:
                    self.bands.append(bands_temp)
                    bands_temp = []
                    k_point+=1
                    # print(k_point, len(self.bands))

        import numpy as np
        numpy_array = np.array(self.bands)
        transpose = numpy_array.T
        self.bands = transpose.tolist()

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
        return self.NumberOfkPoints

    def get_number_of_electrons(self):
        """Returns the total number of electrons in the system"""
        return self.band_index_HOMO

    def get_total_ionic_charge(self):
        """Returns the total ionic charge of the system"""
        return self.band_index_LUMO

    def plot_band_structure(self, **kwargs):
        """Plots the """
        from ebk.BandPlotter import BandPlotter
        import matplotlib.pyplot as plt
        kwargs.update({"arrow_data":[[self.highest_valance[0] , self.highest_valance[1],self.lowest_conduction[0],self.lowest_conduction[1]]]})
        plot = BandPlotter(**kwargs)
        b_gap = self.lowest_conduction[1] - self.highest_valance[1]
        plot.add_to_plot(self.Ef, self.k_dist, self.bands, self.hsp, self.hss, label=f"E$_g$={b_gap: .3f} eV")
        # plt.plot_arrow = [self.highest_valance[0],self.highest_valance[1],self.lowest_conduction[0],self.lowest_conduction[1]]
        # plt.text((0+0.0023)/2+0.004,(self.highest_valance[1]+self.lowest_conduction[1])/2,f"E$_g$={b_gap: .3f} eV")
        plot.plot()

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

    def get_quadrupole_moments(self, file_name = None):
        """
        This part of the code is still under construction and should be moved up when done.


        Calculates the quadrupole for the sytem given.
        Reads in the *.RHO.cube file generated by Denchar.


        The basic structure of the cube file is given here.
        http://paulbourke.net/dataformats/cube/

        Header
        The first two lines of the header are comments, they are generally ignored by parsing packages or used as two default labels.

        The third line has the number of atoms included in the file followed by the position of the origin of the volumetric data.

        The next three lines give the number of voxels along each axis (x, y, z) followed by the axis vector. Note this means the volume need not be aligned with the coordinate axis, indeed it also means it may be sheared although most volumetric packages won't support that. The length of each vector is the length of the side of the voxel thus allowing non cubic volumes. If the sign of the number of voxels in a dimension is positive then the units are Bohr, if negative then Angstroms.

        The last section in the header is one line for each atom consisting of 5 numbers, the first is the atom number, the second is the charge, and the last three are the x,y,z coordinates of the atom center.

        Volumetric dataive the number of voxels along each 
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
        from scipy import constants as c
        print("Calculating the dipoles and quadrupoles...")

        if file_name == None:
            # Here we can default to open the normal output file from Dencahr
            f = open(f"{self.out_file_name}.RHO.cube", "r")
            for line in f:
                self.RHO_file.append(line)
            f.close()
        else:
            read_file_name = file_name
            f = open(f"{read_file_name}", "r")
            for line in f:
                self.RHO_file.append(line)
            f.close()

        origin_0 = []
        number_of_atoms = 0
        a_voxel_counter = 0
        b_voxel_counter = 0
        c_voxel_counter = 0
        full_volume = 0
        a_number_of_voxels = 0
        b_number_of_voxels = 0
        c_number_of_voxels = 0
        a_voxel_length = 0
        b_voxel_length = 0
        c_voxel_length = 0
        a_voxel_vec = 0
        b_voxel_vec = 0
        c_voxel_vec = 0
        original_units = "Angs"
        atom_coordinates_trigger = False
        total_unnormalized_charge = 0
        atoms_info = []# atomic_number, charge, valance_electrons, [x,y,z]
        first_lines_of_file = []
        rrho = 0 # The summation for the r times rho that we need to calculate the centre of charge.
        rrho_elect = 0
        rrho_ionic = 0
        COC = []  # Centre of charge (vector)
        COC_elect = []
        COC_ionic = []

        for i,line in enumerate(self.RHO_file):
            if i == 2:
                first_lines_of_file.append(line)
                parsed = line.split()
                number_of_atoms = int(parsed[0])
                for x in range(1,len(parsed)):
                    origin_0.append(float(parsed[x])) # unit conversions are handled seperalty see below
            if i == 3:
                first_lines_of_file.append(line)
                # print("this is inside i3")
                parsed = line.split()
                a_number_of_voxels = int(parsed[0])
                a_voxel_length = np.sqrt(float(parsed[1])**2+float(parsed[2])**2+float(parsed[3])**2)
                a_voxel_vec = np.array([float(parsed[1]),float(parsed[2]),float(parsed[3])])
                # print(a_voxel_vec)
                # print(parsed)
            if i == 4:
                first_lines_of_file.append(line)
                # print("this is inside i4")
                parsed = line.split()
                b_number_of_voxels = int(parsed[0])
                b_voxel_length = np.sqrt(float(parsed[1])**2+float(parsed[2])**2+float(parsed[3])**2)
                b_voxel_vec = np.array([float(parsed[1]),float(parsed[2]),float(parsed[3])])
                # print(parsed)
            if i == 5:
                first_lines_of_file.append(line)
                # print("this is inside i5")
                parsed = line.split()
                c_number_of_voxels = int(parsed[0])
                c_voxel_length = np.sqrt(float(parsed[1])**2+float(parsed[2])**2+float(parsed[3])**2)
                c_voxel_vec = np.array([float(parsed[1]),float(parsed[2]),float(parsed[3])])
                # print(parsed)
                rho = np.zeros((a_number_of_voxels, b_number_of_voxels, c_number_of_voxels))  # This is the normalized charge density in terms of the number of electrons. It will later be converted into charge in Coulombs.
                rho_ionic = np.zeros((a_number_of_voxels, b_number_of_voxels, c_number_of_voxels))  # This is the normalized charge density in terms of the number of electrons. It will later be converted into charge in Coulombs.

            origin = np.array(origin_0)
            if i == 6: atom_coordinates_trigger = True
            if atom_coordinates_trigger:
                if len(line.split()) != 5:
                    atom_coordinates_trigger = False
            if atom_coordinates_trigger:
                first_lines_of_file.append(line)
                parsed = line.split()
                atomic_number = int(parsed[0])
                charge = float(parsed[1])
                if atomic_number > 10:
                    valance_electrons = atomic_number-10 
                elif atomic_number > 8:
                    valance_electrons = atomic_number-8
                elif atomic_number > 2:
                    valance_electrons = atomic_number-2
                elif atomic_number <= 2: 
                    valance_electrons = atomic_number
                # print(atomic_number, valance_electrons)
                atoms_info.append([atomic_number, charge, valance_electrons, [float(parsed[2]), float(parsed[3]), float(parsed[4])]])
            
            if not atom_coordinates_trigger and i>6:
                parsed = line.split()
                for val in parsed:
                    rho[a_voxel_counter,b_voxel_counter,c_voxel_counter] = float(val)
                    total_unnormalized_charge+= float(val)
                    c_voxel_counter+=1
                    if c_voxel_counter == c_number_of_voxels: 
                        c_voxel_counter = 0
                        b_voxel_counter+=1
                        if b_voxel_counter == b_number_of_voxels:
                            b_voxel_counter=0
                            a_voxel_counter+=1

        # Possible unit conversions are handled here.
        if a_number_of_voxels > 0 and b_number_of_voxels > 0 and c_number_of_voxels > 0:
            original_units = "Bohr"
            # This case means that the units are in Bohr and we have to convert to angstroms
            # COnverting the origin
            for x in range(0,len(origin)):
                    origin[x] = origin[x]*0.529177249
            # Converting the cell vectors into angstroms
            for i,v in enumerate(a_voxel_vec):
                a_voxel_vec[i] = a_voxel_vec[i]*0.529177249
                b_voxel_vec[i] = b_voxel_vec[i]*0.529177249
                c_voxel_vec[i] = c_voxel_vec[i]*0.529177249
            # Converting the atomic coordinates
            for atom in atoms_info:
                # print(atom)
                for i,v in enumerate(a_voxel_vec):
                    atom[3][i] = atom[3][i]*0.529177249
                    # pass
                # print("new",atom)
        d_V = np.cross(a_voxel_vec,b_voxel_vec)
        d_V = np.dot(c_voxel_vec,d_V)

        # Normalizing the charge
        factor = self.number_of_electrons/total_unnormalized_charge
        total_charge = 0
        total_electronic_charge = 0
        for x in range(a_number_of_voxels):
            for y in range(b_number_of_voxels):
                for z in range(c_number_of_voxels):
                    rho[x,y,z] = -rho[x,y,z]*factor
                    total_electronic_charge+=rho[x,y,z]
        rho_elect = rho.copy()

        # Adding the ionic components
        total_ionic_charge = 0
        rrho_ionic_not_binned = 0
        dipole_ionic_not_binned = 0
        voxel_midpoint_vec = 0.5*a_voxel_vec+0.5*b_voxel_vec+0.5*c_voxel_vec
        def get_r_vec(x,y,z):
            # r_vec = ia*a_voxel_vec+ib*b_voxel_vec+ic*c_voxel_vec  + voxel_midpoint_vec
            r_vec = ia*a_voxel_vec+ib*b_voxel_vec+ic*c_voxel_vec + origin
            return r_vec


        prog = progress_bar(len(atoms_info)*a_number_of_voxels*b_number_of_voxels*c_number_of_voxels, descriptor="Loading ions onto grid")
        for i_atom, atom in enumerate(atoms_info):
            # Have to remember that the coordinates are now changed (origin has changed) so we have to use the cube file coordinates
            import math
            coordinates = atom[3]
            # r_atom = np.array([coordinates[0]-origin[0], coordinates[1]-origin[1], coordinates[2]-origin[2]])
            r_atom = np.array([coordinates[0], coordinates[1], coordinates[2]])
            r_voxel = a_voxel_vec+b_voxel_vec+c_voxel_vec

            # # Method 1
            # x_index = math.floor((r_atom[0])/r_voxel[0])
            # y_index = math.floor((r_atom[1])/r_voxel[1])
            # z_index = math.floor((r_atom[2])/r_voxel[2])
            # charge = atom[2]
            # total_ionic_charge += charge
            # rho[x_index, y_index, z_index] += charge
            # rho_ionic[x_index, y_index, z_index] += charge
            # rrho_ionic_not_binned += r_atom*charge
            # dipole_ionic_not_binned += charge*r_atom
            # # prog.get_progress(i_atom)

            # Method 2 takes longer more accurate This is for safe keeping
            for ia in range(a_number_of_voxels):
                for ib in range(b_number_of_voxels):
                    for ic in range(c_number_of_voxels):
                        prog.get_progress(i_atom*a_number_of_voxels*b_number_of_voxels*c_number_of_voxels+(ia)*(b_number_of_voxels*c_number_of_voxels)+(ib)*(c_number_of_voxels)+ic)
                        r_vec = get_r_vec(ia, ib, ic)
                        d_vec = abs(r_atom - r_vec)
                        if d_vec[0] <= voxel_midpoint_vec[0] and d_vec[1] <= voxel_midpoint_vec[1] and d_vec[2] <= voxel_midpoint_vec[2]:
                            charge = atom[2]
                            total_ionic_charge += charge
                            rho[ia, ib, ic] += charge
                            rho_ionic[ia, ib, ic] += charge
                            rrho_ionic_not_binned += r_atom*charge
                            dipole_ionic_not_binned += charge*r_atom
                            break

        print("")
        Qxx = 0
        Qxy = 0
        Qxz = 0
        Qyx = 0
        Qyy = 0
        Qyz = 0
        Qzx = 0
        Qzy = 0
        Qzz = 0
        dipole = 0
        dipole_ionic = 0
        dipole_elect = 0
        total_charge = 0

        # This set of for loops does preliminary calculations that would be required.
        prog = progress_bar(a_number_of_voxels*b_number_of_voxels*c_number_of_voxels, descriptor="Preliminary calculations")
        for ia in range(a_number_of_voxels):
            for ib in range(b_number_of_voxels):
                for ic in range(c_number_of_voxels):
                    prog.get_progress((ia)*(b_number_of_voxels*c_number_of_voxels)+(ib)*(c_number_of_voxels)+ic)
                    """Methana podi indeces proshnayak thiyeanwa. mokenda iterate karanna one i+1 da nattam i da kiyala"""
                    r_vec = get_r_vec(ia, ib, ic)
                    r2 = np.dot(r_vec,r_vec)
                    r = np.sqrt(r2)
                    rrho += rho[ia, ib, ic]*r_vec
                    rrho_ionic += rho_ionic[ia, ib, ic]*r_vec
                    rrho_elect += rho_elect[ia, ib, ic]*r_vec
                    total_charge+=rho[ia, ib, ic]
                    full_volume+=d_V

                    # if ia == 0 and ib == 0 and ic == 0:
                    #     # print("Corner_start", r, r2, r_vec, ia, ib, ic)
                    #     print("Corner_start", r*0.529177249, r2*0.529177249*0.529177249, r_vec*0.529177249, ia*0.529177249, ib*0.529177249, ic*0.529177249)
                    # if ia == 12 and ib == 12 and ic == 12:
                    #     # print("Corner_mid", r, r2, r_vec, ia, ib, ic)
                    #     print("Corner_mid", r*0.529177249, r2*0.529177249*0.529177249, r_vec*0.529177249, ia*0.529177249, ib*0.529177249, ic*0.529177249)
                    # if ia == a_number_of_voxels-1 and ib == b_number_of_voxels-1 and ic == c_number_of_voxels-1:
                    #     print("Corner_end", r*0.529177249, r2*0.529177249*0.529177249, r_vec*0.529177249, ia*0.529177249, ib*0.529177249, ic*0.529177249)

        # Centre of Charge calculations
        COC = rrho/(total_charge)
        COC_elect = rrho_elect/total_electronic_charge
        COC_ionic = rrho_ionic/total_ionic_charge
        COC_ionic_not_binned = rrho_ionic_not_binned/total_ionic_charge
        COC = COC_ionic_not_binned # This should be set so that it can change.

        prog = progress_bar(a_number_of_voxels*b_number_of_voxels*c_number_of_voxels, descriptor="Analyzing for moments")
        for ia in range(a_number_of_voxels):
            for ib in range(b_number_of_voxels):
                for ic in range(c_number_of_voxels):
                    prog.get_progress((ia)*(b_number_of_voxels*c_number_of_voxels)+(ib)*(c_number_of_voxels)+ic)
                    """Methana podi indeces prashnayak thiyeanwa. mokenda iterate karanna one i+1 da nattam i da kiyala"""
                    r_vec = get_r_vec(ia, ib, ic)
                    r_vec = r_vec - COC
                    r2 = np.dot(r_vec,r_vec)
                    r = np.sqrt(r2)
                    dipole+=rho[ia, ib, ic]*r_vec
                    dipole_ionic+=rho_ionic[ia, ib, ic]*r_vec
                    dipole_elect+=rho_elect[ia, ib, ic]*r_vec
                    Qxx+= rho[ia, ib, ic]*(3*(r_vec[0])*(r_vec[0]) - r2)
                    Qyy+= rho[ia, ib, ic]*(3*(r_vec[1])*(r_vec[1]) - r2)
                    Qzz+= rho[ia, ib, ic]*(3*(r_vec[2])*(r_vec[2]) - r2)
                    # Off axis elements
                    Qxy+= rho[ia, ib, ic]*(3*(r_vec[0])*(r_vec[1]))
                    Qxz+= rho[ia, ib, ic]*(3*(r_vec[0])*(r_vec[2]))
                    Qyx+= rho[ia, ib, ic]*(3*(r_vec[1])*(r_vec[0]))
                    Qyz+= rho[ia, ib, ic]*(3*(r_vec[1])*(r_vec[2]))
                    Qzx+= rho[ia, ib, ic]*(3*(r_vec[2])*(r_vec[0]))
                    Qzy+= rho[ia, ib, ic]*(3*(r_vec[2])*(r_vec[1]))

        # Calculations happen here
        # Printing values out for now. Should make proper get methods.
        print("Origin", origin)
        print("numer of voxels", a_number_of_voxels, b_number_of_voxels, c_number_of_voxels)
        print("Full volume", full_volume)
        print("Half cell vector", voxel_midpoint_vec)
        print("Voxel cell vector", r_voxel)
        print("elementary volume", d_V)
        print("elementary charge scipy", c.elementary_charge)
        print("Volume we should get", 25*1.023602*25*1.023602*25*1.653511 )
        print("check for debye", c.elementary_charge*10**(-10)/(3.336*10**(-30)))
        print("Denchar Volume", 13*13*21 )
        print("dipole", dipole*c.elementary_charge*10**(-10)/(3.336*10**(-30)), "Debye")
        print("dipole elec", dipole_elect*c.elementary_charge*10**(-10)/(3.336*10**(-30)), "Debye")
        print("dipole ionic", dipole_ionic*c.elementary_charge*10**(-10)/(3.336*10**(-30)), "Debye")
        print("dipole non binned", dipole_ionic_not_binned*c.elementary_charge*10**(-10)/(3.336*10**(-30)), "Debye")
        print("total_electrons", self.number_of_electrons)
        print("normalized electronic charge", total_electronic_charge)
        print("total ionic charge", total_ionic_charge)
        print("Total Charge", total_charge)
        print("Qxx", Qxx*c.elementary_charge*10**(-10)/(3.336*10**(-30)), "Debye")
        print("Qyy", Qyy*c.elementary_charge*10**(-10)/(3.336*10**(-30)), "Debye")
        print("Qzz", Qzz*c.elementary_charge*10**(-10)/(3.336*10**(-30)), "Debye")
        # print("rrho rhoelectronic and rrho ionic", rrho, rho_elect, rho_ionic)
        print("shapes", np.shape(rrho), np.shape(rrho_elect), np.shape(rrho_ionic))
        print("Centre of Charge (electronic, ionic, ionic_binned)", COC_elect, COC_ionic, COC_ionic_not_binned)
        print("Centre of Charge (ionic - elec)", COC_ionic - COC_elect)

        # Writing out the new data in a new file
        with open(f"{self.out_file_name}.electronicandionic.cube", "w+") as file:
            file.write("Produced by Eranjan\n")
            file.write(f"For total charge for system: {self.out_file_name}\n")
            for line in first_lines_of_file:
                file.write(line)
            for x in range(a_number_of_voxels):
                for y in range(b_number_of_voxels):
                    for z in range(c_number_of_voxels):            
                        file.write(f"{rho[x][y][z]:1.5E} ")
                        if (z % 6 == 5):
                         file.write("\n")
                    file.write("\n")
    

        # Writing out the new data in a new file
        with open(f"{self.out_file_name}.ionic.cube", "w+") as file:
            file.write("Produced by Eranjan\n")
            file.write(f"For ionic charge for system: {self.out_file_name}\n")
            for line in first_lines_of_file:
                file.write(line)
            for x in range(a_number_of_voxels):
                for y in range(b_number_of_voxels):
                    for z in range(c_number_of_voxels):            
                        file.write(f"{rho_ionic[x][y][z]:1.5E} ")
                        if (z % 6 == 5):
                         file.write("\n")
                    file.write("\n")
    
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
