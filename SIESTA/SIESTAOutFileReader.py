"""
This file reads the the file out_file_name.out and extracts/calculates data/values from it
"""
from ebk import progress_bar
from numpy import char
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
        trigger_read_in_cell = False

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

            if "siesta: Electric dipole (Debye) =" in line:
                x = line.strip("siesta: Electric dipole (Debye) =")
                x = x.split()
                self.dipole_siesta = [float(x[0]), float(x[1]), float(x[2])]

            # Here we are reading in the Total energies
            if "siesta:         Total =" in line:
                x = line.strip("siesta:         Total =")
                x = x.split()
                self.E_total = float(x[0])      

            # Here we are reading in the Total energies
            if "siesta: Automatic unit cell vectors (Ang):" in line:
                trigger_read_in_cell = True
                self.initial_cell_vectors = []
                continue
            if trigger_read_in_cell:
                x = line.strip("siesta:")
                x = x.split()
                x = [float(i) for i in x]
                self.initial_cell_vectors.append(x)
                if len(self.initial_cell_vectors) == 3: trigger_read_in_cell = False

            # Here we are reading in the number of atoms and orbitals and projectors
            if "initatomlists: Number of atoms, orbitals, and projectors:" in line:
                x = line.strip("initatomlists: Number of atoms, orbitals, and projectors:")
                x = x.split()
                self.N_atoms = int(x[0])
                self.N_orbitals = int(x[1])
                self.N_projectors = int(x[2])

            # Here we are reading in the Iternal realspace grid
            if "InitMesh: (bp) = " in line:
                x = line.strip("InitMesh: (bp) = ")
                x = x.split()
                self.InitMesh_bp = [int(x[0]), int(x[2]), int(x[4])]

            # Here we are reading in the External realspace grid
            if "ExtMesh (bp) on 0 = " in line:
                x = line.strip("ExtMesh (bp) on 0 = ")
                x = x.split()
                self.ExtMesh_bp = [int(x[0]), int(x[2]), int(x[4])]

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
        """Calculates the index of the HOMO level. The indeces start from 1"""
        return self.band_index_HOMO

    def get_HOMO_energy_for_gammapoint_calculation(self):
        """This only works for the gamma point since otherwise we would have to to choose the right k point at gamma"""
        return self.eigen_values_by_k_point[1][self.get_HOMO_band_index() - 1]

    def get_LUMO_energy_for_gammapoint_calculation(self):
        """This only works for the gamma point since otherwise we would have to to choose the right k point at gamma"""
        return self.eigen_values_by_k_point[1][self.get_LUMO_band_index() - 1]

    def get_LUMO_band_index(self):
        """Calculates the index of the LUMO level. The indeces start from 1"""
        return self.band_index_LUMO

    def get_number_of_electrons(self):
        """Returns the total number of electrons in the system"""
        return self.number_of_electrons

    def get_total_ionic_charge(self):
        """Returns the total ionic charge of the system"""
        return self.total_ionic_charge

    def get_initial_cell_vectors(self):
        """Returns the un relaxed initial cell vectors"""
        return self.initial_cell_vectors

    def plot_band_structure(self, **kwargs):
        """Plots the """
        from ebk.BandPlotter import BandPlotter
        # import matplotlib.pyplot as plt
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

    def get_quadrupole_moments(self, file_name = None, out_put_file_name = "", cell = [[0., 0., 0.,],[0., 0., 0.,],[0., 0., 0.,]]):
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
        print("\nCalculating the dipoles and quadrupoles...")

        if file_name == None:
            # Here we can default to open the normal output file from Dencahr
            read_file_name = f"{self.out_file_name}.RHO.cube"
        else:
            read_file_name = file_name
        print(f"Reading in file {read_file_name}")   # This is important so that the user knows which file is being read
        f = open(f"{read_file_name}", "r")
        for line in f:
            self.RHO_file.append(line)
        f.close()

        number_of_atoms = 0
        full_volume = 0
        a_voxel_counter = 0
        b_voxel_counter = 0
        c_voxel_counter = 0
        original_units = "Angs"
        atom_coordinates_trigger = False
        total_unnormalized_charge = 0
        atoms_info = []# atomic_number, charge, valance_electrons, [x,y,z]
        first_lines_of_file = []

        for i,line in enumerate(self.RHO_file):
            if i == 2:
                first_lines_of_file.append(line)
                parsed = line.split()
                number_of_atoms = int(parsed[0])
                atoms_left_to_add = number_of_atoms
                origin = np.array([float(parsed[1]),float(parsed[2]),float(parsed[3])]) # unit conversions are handled seperalty see below
            if i == 3:
                first_lines_of_file.append(line)
                parsed = line.split()
                a_number_of_voxels = int(parsed[0])
                a_voxel_vec = np.array([float(parsed[1]),float(parsed[2]),float(parsed[3])])
            if i == 4:
                first_lines_of_file.append(line)
                parsed = line.split()
                b_number_of_voxels = int(parsed[0])
                b_voxel_vec = np.array([float(parsed[1]),float(parsed[2]),float(parsed[3])])
            if i == 5:
                first_lines_of_file.append(line)
                parsed = line.split()
                c_number_of_voxels = int(parsed[0])
                c_voxel_vec = np.array([float(parsed[1]),float(parsed[2]),float(parsed[3])])
                rho = np.zeros((a_number_of_voxels, b_number_of_voxels, c_number_of_voxels))  # This is the normalized charge density in terms of the number of electrons. It will later be converted into charge in Coulombs.
            
            # Here we load in the atoms from the cube file too
            if i == 6: atom_coordinates_trigger = True
            if atom_coordinates_trigger:
                if atoms_left_to_add == 0:
                    atom_coordinates_trigger = False
            if atom_coordinates_trigger:
                first_lines_of_file.append(line)
                parsed = line.split()
                atomic_number = int(parsed[0])
                charge = float(parsed[1])            # This is not used anywhere in the code and is not needed for our cases.
                atoms_left_to_add-=1
                if atomic_number == 1000:            # For testing code
                    valance_electrons = 0.058874798240219241 
                elif atomic_number == 1:               # H
                    valance_electrons = 1 
                elif atomic_number == 6:               # O
                    valance_electrons = 4 
                elif atomic_number == 7:               # N
                    valance_electrons = 5 
                elif atomic_number == 8:
                    valance_electrons = 6
                elif atomic_number == 9:
                    valance_electrons = 7 
                elif atomic_number == 16:              # S
                    valance_electrons = 6 
                elif atomic_number == 17:              # Cl
                    valance_electrons = 7 
                elif atomic_number == 26:
                    valance_electrons = 7 
                elif atomic_number == 50:              # Sn
                    valance_electrons = 4 
                else:
                    print("get_quadrupole_moments: Atoms not recognized")
                    return
                atoms_info.append([atomic_number, charge, valance_electrons, [float(parsed[2]), float(parsed[3]), float(parsed[4])]])
          
            # Here we load the volumetric data
            if not atom_coordinates_trigger and i>5:
                parsed = line.split()
                for val in parsed:
                    val_read = float(val)
                    rho[a_voxel_counter, b_voxel_counter, c_voxel_counter] = val_read
                    total_unnormalized_charge+= val_read
                    c_voxel_counter+=1
                    if c_voxel_counter == c_number_of_voxels: 
                        c_voxel_counter = 0
                        b_voxel_counter+=1
                        if b_voxel_counter == b_number_of_voxels:
                            b_voxel_counter=0
                            a_voxel_counter+=1
        # All data is loaded by this point

        # Possible unit conversions are handled here.
        if a_number_of_voxels > 0 and b_number_of_voxels > 0 and c_number_of_voxels > 0:
            print("Units are detected to be Bohr and are converted into Angstroms")
            # This case means that the units are in Bohr and we have to convert to angstroms
            original_units = "Bohr"
            # Converting the origin and cell vectors into angstroms
            origin = origin*0.529177249
            a_voxel_vec = a_voxel_vec*0.529177249
            b_voxel_vec = b_voxel_vec*0.529177249
            c_voxel_vec = c_voxel_vec*0.529177249
            # Converting the atomic coordinates
            for atom in atoms_info:
                for i,v in enumerate(a_voxel_vec):
                    atom[3][i] = atom[3][i]*0.529177249

        # Here some elemetary volume and distances are computed
        r_voxel = a_voxel_vec+b_voxel_vec+c_voxel_vec
        d_V = np.dot(c_voxel_vec, np.cross(a_voxel_vec,b_voxel_vec))

        # Normalizing the charge
        if total_unnormalized_charge !=0:
            factor = self.number_of_electrons/total_unnormalized_charge  # To get around the devide by zero error for the zero charge case
            # factor = 1
        else: factor = 0   
        total_electronic_charge = 0
        for x in range(a_number_of_voxels):
            for y in range(b_number_of_voxels):
                for z in range(c_number_of_voxels):
                    rho[x,y,z] = -rho[x,y,z]*factor
                    total_electronic_charge+=rho[x,y,z]

        # Adding the ionic components and the dipole part for the the non binned ionic part
        total_ionic_charge = 0
        voxel_midpoint_vec = 0.5*a_voxel_vec+0.5*b_voxel_vec+0.5*c_voxel_vec
        def get_r_vec(x,y,z):
            r_vec = x*a_voxel_vec+y*b_voxel_vec+z*c_voxel_vec + origin
            return r_vec

        Qxx = 0
        Qxy = 0
        Qxz = 0
        Qyx = 0
        Qyy = 0
        Qyz = 0
        Qzx = 0
        Qzy = 0
        Qzz = 0
        Qxx_non_traceless = 0
        Qxy_non_traceless = 0
        Qxz_non_traceless = 0
        Qyx_non_traceless = 0
        Qyy_non_traceless = 0
        Qyz_non_traceless = 0
        Qzx_non_traceless = 0
        Qzy_non_traceless = 0
        Qzz_non_traceless = 0
        total_charge = 0
        dipole = 0

        prog = progress_bar(a_number_of_voxels*b_number_of_voxels*c_number_of_voxels, descriptor="Analyzing for moments")
        for ia in range(a_number_of_voxels):
            for ib in range(b_number_of_voxels):
                for ic in range(c_number_of_voxels):
                    prog.get_progress((ia)*(b_number_of_voxels*c_number_of_voxels)+(ib)*(c_number_of_voxels)+ic)
                    """Methana podi indeces prashnayak thiyeanwa. mokenda iterate karanna one i+1 da nattam i da kiyala"""
                    r_vec0 = get_r_vec(ia, ib, ic)
                    r_vec = r_vec0
                    r2 = np.dot(r_vec,r_vec)
                    r = np.sqrt(r2)
                    total_charge+=rho[ia, ib, ic]
                    full_volume+=d_V
                    dipole+=rho[ia, ib, ic]*r_vec

                    # For the traceless part
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
                    # Now for the non-traceless part
                    Qxx_non_traceless+= rho[ia, ib, ic]*(r_vec[0])*(r_vec[0])
                    Qyy_non_traceless+= rho[ia, ib, ic]*(r_vec[1])*(r_vec[1])
                    Qzz_non_traceless+= rho[ia, ib, ic]*(r_vec[2])*(r_vec[2])
                    # Off axis elements
                    Qxy_non_traceless+= rho[ia, ib, ic]*(r_vec[0])*(r_vec[1])
                    Qxz_non_traceless+= rho[ia, ib, ic]*(r_vec[0])*(r_vec[2])
                    Qyx_non_traceless+= rho[ia, ib, ic]*(r_vec[1])*(r_vec[0])
                    Qyz_non_traceless+= rho[ia, ib, ic]*(r_vec[1])*(r_vec[2])
                    Qzx_non_traceless+= rho[ia, ib, ic]*(r_vec[2])*(r_vec[0])
                    Qzy_non_traceless+= rho[ia, ib, ic]*(r_vec[2])*(r_vec[1])

        # This part adds in the ionic non binned part to the quadrupoles and the dipoles
        for i_atom, atom in enumerate(atoms_info):
            coordinates = atom[3]
            charge = atom[2]   # This is from the valance electron item in the list
            # print("charge",charge)
            r_atom = np.array([coordinates[0], coordinates[1], coordinates[2]])
            r_vec = r_atom
            r2 = np.dot(r_vec,r_vec)
            dipole += charge*r_vec
            Qxx+= charge*(3*(r_vec[0])*(r_vec[0]) - r2)
            Qyy+= charge*(3*(r_vec[1])*(r_vec[1]) - r2)
            Qzz+= charge*(3*(r_vec[2])*(r_vec[2]) - r2)
            # Off axis elements
            Qxy+= charge*(3*(r_vec[0])*(r_vec[1]))
            Qxz+= charge*(3*(r_vec[0])*(r_vec[2]))
            Qyx+= charge*(3*(r_vec[1])*(r_vec[0]))
            Qyz+= charge*(3*(r_vec[1])*(r_vec[2]))
            Qzx+= charge*(3*(r_vec[2])*(r_vec[0]))
            Qzy+= charge*(3*(r_vec[2])*(r_vec[1]))
            # Now for the non-tracelss part
            Qxx_non_traceless+= charge*(r_vec[0])*(r_vec[0])
            Qyy_non_traceless+= charge*(r_vec[1])*(r_vec[1])
            Qzz_non_traceless+= charge*(r_vec[2])*(r_vec[2])
            # print("Qzz_atom",Qzz_non_traceless, charge, r_vec[2], charge*(r_vec[2])*(r_vec[2]))
            # Off axis elements
            Qxy_non_traceless+= charge*(r_vec[0])*(r_vec[1])
            Qxz_non_traceless+= charge*(r_vec[0])*(r_vec[2])
            Qyx_non_traceless+= charge*(r_vec[1])*(r_vec[0])
            Qyz_non_traceless+= charge*(r_vec[1])*(r_vec[2])
            Qzx_non_traceless+= charge*(r_vec[2])*(r_vec[0])
            Qzy_non_traceless+= charge*(r_vec[2])*(r_vec[1])
            total_ionic_charge+=charge
            total_charge += charge

        # Some units
        unit_factor_Debye = 1/0.2081943 #(same as c.elementary_charge*10**(-10)/(3.33564*10**(-30)))
        # unit_factor_Debye = 1*10**(10)/(3.33564*10**(-30))
        unit_factor_esuA = unit_factor_Debye*10**(-10)
        unit_factor_Cmm = c.elementary_charge*10**(-20)
        conversions = {"Debye.Angs":unit_factor_Debye, "Cm^2":unit_factor_Cmm,"esu.Angs":unit_factor_esuA, "in q.r^2 numofelectrons.angs^2":1}

        # Here we can write to output file and then delete the above if necessary for furture 
        if out_put_file_name == "": out_put_file_name = self.out_file_name
        print(f"Writing the values to file {out_put_file_name}.electrostatics.out")
        summary_file = open(f"{out_put_file_name}.electrostatics.out", "w")
        summary_file.write(f"Summary of electrostatics calculations for file {self.out_file_name}\n\n")
        summary_file.write(f"numer of voxels                            {a_number_of_voxels}, {b_number_of_voxels}, {c_number_of_voxels}\n")
        summary_file.write(f"Voxel vectors                              {a_voxel_vec}, {b_voxel_vec}, {c_voxel_vec}\n")
        summary_file.write(f"Origin                                     {origin}\n")
        summary_file.write(f"Original Units                             {original_units}\n")
        summary_file.write(f"Full volume                                {full_volume:<5f}\n")
        summary_file.write(f"Half cell vector                           {voxel_midpoint_vec}\n")
        summary_file.write(f"Voxel cell vector                          {r_voxel}\n")
        summary_file.write(f"elementary volume                          {d_V:<5f}\n")
        summary_file.write(f"total_electrons (from siesta)              {self.number_of_electrons}\n")
        summary_file.write(f"normalized total electronic charge         {total_electronic_charge}\n")
        summary_file.write(f"total ionic charge                         {total_ionic_charge}\n")
        summary_file.write(f"Normalization factor for electron density  {factor:<5f}\n")
        summary_file.write(f"Total Unnormalized charge from cube file   {total_unnormalized_charge:<5f}\n")
        summary_file.write(f"Total Charge                               {total_charge}\n")
        summary_file.write(f"dipole                                     {dipole*unit_factor_Debye} Debye\n")
        summary_file.write(f"dipole                                     {dipole} in q.r numofelectrons.angs\n")
        for i,(k,v)in enumerate(conversions.items()):
            summary_file.write(f"Qxx                                        {Qxx*v:<5f} {k}\n")
            summary_file.write(f"Qyy                                        {Qyy*v:<5f} {k}\n")
            summary_file.write(f"Qzz                                        {Qzz*v:<5f} {k}\n")
            summary_file.write(f"Qxy                                        {Qxy*v:<5f} {k}\n")
            summary_file.write(f"Qxz                                        {Qxz*v:<5f} {k}\n")
            summary_file.write(f"Qyx                                        {Qyx*v:<5f} {k}\n")
            summary_file.write(f"Qyz                                        {Qyz*v:<5f} {k}\n")
            summary_file.write(f"Qzx                                        {Qzx*v:<5f} {k}\n")
            summary_file.write(f"Qzy                                        {Qzy*v:<5f} {k}\n")
            summary_file.write(f"Qxx  (in the non-traceless from)           {Qxx_non_traceless*v:<5f} {k}\n")
            summary_file.write(f"Qyy  (in the non-traceless from)           {Qyy_non_traceless*v:<5f} {k}\n")
            summary_file.write(f"Qzz  (in the non-traceless from)           {Qzz_non_traceless*v:<5f} {k}\n")
            summary_file.write(f"Qxy  (in the non-traceless from)           {Qxy_non_traceless*v:<5f} {k}\n")
            summary_file.write(f"Qxz  (in the non-traceless from)           {Qxz_non_traceless*v:<5f} {k}\n")
            summary_file.write(f"Qyx  (in the non-traceless from)           {Qyx_non_traceless*v:<5f} {k}\n")
            summary_file.write(f"Qyz  (in the non-traceless from)           {Qyz_non_traceless*v:<5f} {k}\n")
            summary_file.write(f"Qzx  (in the non-traceless from)           {Qzx_non_traceless*v:<5f} {k}\n")
            summary_file.write(f"Qzy  (in the non-traceless from)           {Qzy_non_traceless*v:<5f} {k}\n")
        # summary_file.write(f"C.O.C (if no ions then binned electronic)  {COC}\n")
        summary_file.close()

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
    
        return_dict = {}
        return_dict.update({"Qzz":Qzz*unit_factor_Debye})
        print("Electrostatics: Done\n")

        return return_dict        

    def load_quadrupole_moments(self, file_name = None, cell = [[0., 0., 0.,],[0., 0., 0.,],[0., 0., 0.,]]):
        """
        Code reads in the Filename.electrostatics.out file generated by the get_quadupole_moments

        dipoles are read in with units Debye
        quadrupoles are read in with units Debye.Angs

        """
        if file_name == None:
            # Here we can default to open the normal output file from Dencahr
            read_file_name = f"{self.out_file_name}.electrostatics.out"
        else:
            read_file_name = file_name
        print(f"Reading in file {read_file_name}")
        f = open(f"{read_file_name}", "r")
        for line in f:
            if "dipole" in line and "Debye" in line and "unadjusted" not in line and "binned" not in line:
                line = line.strip("[")
                line = line.strip("]")
                # print(line)
                parsed = line.split()
                if "[" in parsed : parsed.remove("[")
                if "]" in parsed : parsed.remove("]")
                # print(parsed)
                self.dipole = [ float(parsed[-4].strip("[")), float(parsed[-3]), float(parsed[-2].strip("]")) ]
            if "Qxx" in line and "Debye.Angs" in line and "(in the non-traceless from)" not in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Qxx = float(parsed[1])
            if "Qyy" in line and "Debye.Angs" in line and "(in the non-traceless from)" not in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Qyy = float(parsed[1])
            if "Qzz" in line and "Debye.Angs" in line and "(in the non-traceless from)" not in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Qzz = float(parsed[1])
            if "Qxy" in line and "Debye.Angs" in line and "(in the non-traceless from)" not in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Qxy = float(parsed[1])
            if "Qxz" in line and "Debye.Angs" in line and "(in the non-traceless from)" not in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Qxz = float(parsed[1])
            if "Qyx" in line and "Debye.Angs" in line and "(in the non-traceless from)" not in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Qyx = float(parsed[1])
            if "Qyz" in line and "Debye.Angs" in line and "(in the non-traceless from)" not in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Qyz = float(parsed[1])
            if "Qzx" in line and "Debye.Angs" in line and "(in the non-traceless from)" not in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Qzx = float(parsed[1])
            if "Qzy" in line and "Debye.Angs" in line and "(in the non-traceless from)" not in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Qzy = float(parsed[1])

            if "Qxx" in line and "Debye.Angs" in line and "(in the non-traceless from)" in line and "numofelectrons.angs^2" not in line:
                parsed = line.split() 
                # print(parsed)
                self.Qxx_non_traceless = float(parsed[5])
            if "Qyy" in line and "Debye.Angs" in line and "(in the non-traceless from)" in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Qyy_non_traceless = float(parsed[5])
            if "Qzz" in line and "Debye.Angs" in line and "(in the non-traceless from)" in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Qzz_non_traceless = float(parsed[5])
            if "Qxy" in line and "Debye.Angs" in line and "(in the non-traceless from)" in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Qxy_non_traceless = float(parsed[5])
            if "Qxz" in line and "Debye.Angs" in line and "(in the non-traceless from)" in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Qxz_non_traceless = float(parsed[5])
            if "Qyx" in line and "Debye.Angs" in line and "(in the non-traceless from)" in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Qyx_non_traceless = float(parsed[5])
            if "Qyz" in line and "Debye.Angs" in line and "(in the non-traceless from)" in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Qyz_non_traceless = float(parsed[5])
            if "Qzx" in line and "Debye.Angs" in line and "(in the non-traceless from)" in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Qzx_non_traceless = float(parsed[5])
            if "Qzy" in line and "Debye.Angs" in line and "(in the non-traceless from)" in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Qzy_non_traceless = float(parsed[5])

            if "C.O.C (if no ions then binned electronic)" in line:
                parsed = line.split()
                self.COC = [ float(parsed[-3].strip("[")), float(parsed[-2]), float(parsed[-1].strip("]")) ]
        f.close()

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
