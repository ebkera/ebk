"""
This file reads the the file out_file_name.out and extracts/calculates data/values from it
"""
import imp
from ebk import progress_bar
import numpy as np

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
                self.total_ionic_charge_siesta = float(line.split()[3])

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
        return self.total_ionic_charge_siesta

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


    def generate_new_denchar_file(self, file_name=None, out_file_name=None, number_of_points=[None, None,None], box_edges_input=["Mid","Mid","Mid"]):
        """
        Since it is important to readjust cells and also recompute the charge density this script 
        readin current Denchar files and write a new dencchar file based off of the original

        This can be implemented if needed but not implemented now
        runs_mode: The method creates a .sh file to run all the newly created files
                    THis can be done in 2 modes and append mode and an new file 

        This can be implemented if needed but not implemented now
        runs_mode: The method creates a .sh file to run all the newly created files
                    THis can be done in 2 modes and append mode and an new file 


        box_edges: The type should be like this [[xmin,xmax],[ymin,ymax],[zmin,zmax]]
        """

        if file_name == None:
            # Here we can default to open the normal output file from Denchar
            read_file_name = f"{self.out_file_name}.Denchar.fdf"
            print("Reading in default denchar file...")
        else:
            read_file_name = file_name

        init_cell = self.get_initial_cell_vectors()
        box_edges = [[[],[]],[[],[]],[[],[]]]
        if "Mid" in box_edges_input[0]:
            # We take box edges to be the middle of the images 
            box_edges[0][0] = -(init_cell[0][0]/2)
            box_edges[0][1] = (init_cell[0][0]/2)
        else:box_edges[0] = box_edges_input[0]
        if "Mid" in box_edges_input[1]:
            # We take box edges to be the middle of the images 
            box_edges[1][0] = -(init_cell[1][1]/2)
            box_edges[1][1] = (init_cell[1][1]/2)
        else:box_edges[1] = box_edges_input[1]
        if "Mid" in box_edges_input[2]:
            # We take box edges to be the middle of the images 
            box_edges[2][0] = -(init_cell[2][2]/2)
            box_edges[2][1] = (init_cell[2][2]/2)
        else:box_edges[2] = box_edges_input[2]

        if None in number_of_points:
            print(f"generate_denchar_file: number of bins automatically set to InitMesh_bp: {self.InitMesh_bp}")
            number_of_points = self.InitMesh_bp
            
        if out_file_name == None:
            # Here we can default to we default to using a known file name for denchar
            new_denchar_file_name = f"{self.out_file_name}.newDenchar.{number_of_points}.{box_edges}.fdf"
        else:
            new_denchar_file_name = file_name

        bash_file = open(f"{self.out_file_name}.denchar_runs.sh", "a+")
        new_denchar = open(new_denchar_file_name, "w+")
           
        with open(read_file_name, "r") as old_denchar:
            for line in old_denchar:
                if "Denchar.M" not in line and "NumberPoints" not in line:
                    new_denchar.write(line)

        new_denchar.write(f"Denchar.NumberPointsX       {number_of_points[0]} \n")             
        new_denchar.write(f"Denchar.NumberPointsY       {number_of_points[1]} \n")             
        new_denchar.write(f"Denchar.NumberPointsZ       {number_of_points[2]} \n")
            
        new_denchar.write(f"Denchar.MinX    {box_edges[0][0]} Ang\n")    
        new_denchar.write(f"Denchar.MaxX    {box_edges[0][1]} Ang\n")    
        new_denchar.write(f"Denchar.MinY    {box_edges[1][0]} Ang\n")    
        new_denchar.write(f"Denchar.MaxY    {box_edges[1][1]} Ang\n")    
        new_denchar.write(f"Denchar.MinZ    {box_edges[2][0]} Ang\n")    
        new_denchar.write(f"Denchar.MaxZ    {box_edges[2][1]} Ang\n") 

        new_denchar.close()
        bash_file.write(f"denchar < 'CO.newDenchar.{number_of_points}.{box_edges}.fdf' | tee 'CO.denchar.{number_of_points}.{box_edges}.out'\n")
        bash_file.write(f"cp CO.RHO.cube 'CO.RHO.{number_of_points}.{box_edges}.cube'\n")
        bash_file.close()

    def read_in_rho_cube_file(self, file_name = None):
        """
        Reads in desired *.RHO.cube file generated by Denchar and loads it into a grid (matrix)
        It also calculates some of the basic calculations that are essential.

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
        a_voxel_counter = 0
        b_voxel_counter = 0
        c_voxel_counter = 0
        atom_coordinates_trigger = False
        self.atoms_info = []# atomic_number, charge, valance_electrons, [x,y,z]
        self.first_lines_of_file = []

        for i,line in enumerate(self.RHO_file):
            if i == 2:
                self.first_lines_of_file.append(line)
                parsed = line.split()
                number_of_atoms = int(parsed[0])
                atoms_left_to_add = number_of_atoms
                self.origin = np.array([float(parsed[1]),float(parsed[2]),float(parsed[3])]) # unit conversions are handled seperalty see below
            if i == 3:
                self.first_lines_of_file.append(line)
                parsed = line.split()
                self.a_number_of_voxels = int(parsed[0])
                self.a_voxel_vec = np.array([float(parsed[1]),float(parsed[2]),float(parsed[3])])
            if i == 4:
                self.first_lines_of_file.append(line)
                parsed = line.split()
                self.b_number_of_voxels = int(parsed[0])
                self.b_voxel_vec = np.array([float(parsed[1]),float(parsed[2]),float(parsed[3])])
            if i == 5:
                self.first_lines_of_file.append(line)
                parsed = line.split()
                self.c_number_of_voxels = int(parsed[0])
                self.c_voxel_vec = np.array([float(parsed[1]),float(parsed[2]),float(parsed[3])])
                self.volumetric_data = np.zeros((self.a_number_of_voxels, self.b_number_of_voxels, self.c_number_of_voxels))  # This is the normalized charge density in terms of the number of electrons. It will later be converted into charge in Coulombs.
            
            # Here we load in the atoms from the cube file too
            if i == 6: atom_coordinates_trigger = True
            if atom_coordinates_trigger:
                if atoms_left_to_add == 0:
                    atom_coordinates_trigger = False
            if atom_coordinates_trigger:
                self.first_lines_of_file.append(line)
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
                    print("ERROR: get_quadrupole_moments: Atoms not recognized")
                    return
                self.atoms_info.append([atomic_number, charge, valance_electrons, [float(parsed[2]), float(parsed[3]), float(parsed[4])]])
          
            # Here we load the volumetric data
            if not atom_coordinates_trigger and i>5:
                parsed = line.split()
                for val in parsed:
                    val_read = float(val)
                    self.volumetric_data[a_voxel_counter, b_voxel_counter, c_voxel_counter] = val_read
                    c_voxel_counter+=1
                    if c_voxel_counter == self.c_number_of_voxels: 
                        c_voxel_counter = 0
                        b_voxel_counter+=1
                        if b_voxel_counter == self.b_number_of_voxels:
                            b_voxel_counter=0
                            a_voxel_counter+=1
        # All data is loaded by this point

        # Possible unit conversions are handled here.   
        if self.a_number_of_voxels > 0 and self.b_number_of_voxels > 0 and self.c_number_of_voxels > 0:
            print("Units are detected to be Bohr and are converted into Angstroms")
            # This case means that the units are in Bohr and we have to convert to angstroms
            self.original_units = "Bohr"
            # Converting the origin and cell vectors into angstroms
            self.origin = self.origin*0.529177249
            self.a_voxel_vec = self.a_voxel_vec*0.529177249
            self.b_voxel_vec = self.b_voxel_vec*0.529177249
            self.c_voxel_vec = self.c_voxel_vec*0.529177249
            # Converting the atomic coordinates
            for atom in self.atoms_info:
                for i,v in enumerate(self.a_voxel_vec):
                    atom[3][i] = atom[3][i]*0.529177249
        else: self.original_units = "Angs" 

    def normalize_electronic_charge_density(self):
        """Here this is exclusively for initializations of the inputs read from the read_in_rho_cube_file() method"""
        # Normalizing the charge
        self.rho = np.zeros((self.a_number_of_voxels, self.b_number_of_voxels, self.c_number_of_voxels))  # This is the normalized charge density in terms of the number of electrons. It will later be converted into charge in Coulombs.
        self.total_unnormalized_charge = 0
        self.total_normalized_electronic_charge = 0

        print("Normalizing Charge density please wait...  ", end="", flush=True)
        self.get_volumes_from_cube_header()     # To get dV

        # Here we integrate the total charge w.r.t. the volume.
        for x in range(self.a_number_of_voxels):
            for y in range(self.b_number_of_voxels):
                for z in range(self.c_number_of_voxels):
                    self.total_unnormalized_charge += self.volumetric_data[x,y,z]*self.d_V

        if self.total_unnormalized_charge !=0:
            self.charge_normalization_factor = self.number_of_electrons/self.total_unnormalized_charge  # To get around the devide by zero error for the zero charge case
        else: self.charge_normalization_factor = 0   

        # Here we integrate the total charge w.r.t. the volume.
        for x in range(self.a_number_of_voxels):
            for y in range(self.b_number_of_voxels):
                for z in range(self.c_number_of_voxels):
                    self.rho[x,y,z] = -self.volumetric_data[x,y,z]*self.charge_normalization_factor
                    self.total_normalized_electronic_charge += self.rho[x,y,z]*self.charge_normalization_factor*self.d_V

        print("Done")
 
    def get_volumes_from_cube_header(self):
        """"""
        # Here some initializations and calculations
        self.d_V = np.dot(self.c_voxel_vec, np.cross(self.a_voxel_vec,self.b_voxel_vec))
        self.volume_from_cube_header=np.dot(self.c_voxel_vec*self.c_number_of_voxels, np.cross(self.a_voxel_vec*self.a_number_of_voxels,self.b_voxel_vec*self.b_number_of_voxels)) 

    def write_to_quadrupole_outputfiles(self, out_put_file_name):
        """
        This method writes to outputfiles
        """
        from scipy import constants as c

        def get_direction(number):
            if number == 0:return "x"
            if number == 1:return "y"
            if number == 2:return "z"

        # Writing to the output file
        # Some units and unit conversions
        unit_factor_Debye = 1/0.2081943 #(same as c.elementary_charge*10**(-10)/(3.33564*10**(-30)))
        # unit_factor_Debye = 1*10**(10)/(3.33564*10**(-30))
        unit_factor_esuA = unit_factor_Debye*10**(-10)
        unit_factor_Cmm = c.elementary_charge*10**(-20)
        conversions = {"Debye.Angs":unit_factor_Debye, "Cm^2":unit_factor_Cmm,"esu.Angs":unit_factor_esuA, "in q.r^2 numofelectrons.angs^2":1}

        # Here we can write to output file 
        if out_put_file_name == "": out_put_file_name = self.out_file_name
        print(f"Writing the values to file {out_put_file_name}.electrostatics.out")
        summary_file = open(f"{out_put_file_name}.electrostatics.out", "w")
        summary_file.write(f"Summary of electrostatics calculations for file {self.out_file_name}\n\n")
        summary_file.write(f"numer of voxels                            {self.a_number_of_voxels}, {self.b_number_of_voxels}, {self.c_number_of_voxels}\n")
        summary_file.write(f"Voxel vectors                              {self.a_voxel_vec}, {self.b_voxel_vec}, {self.c_voxel_vec}\n")
        summary_file.write(f"Origin                                     {self.origin}\n")
        summary_file.write(f"Original Units                             {self.original_units}\n")
        summary_file.write(f"Full volume (from cube header)             {self.volume_from_cube_header:<5f}\n")
        if hasattr(self,"summed_volume"):summary_file.write(f"Full volume (summed)                       {self.summed_volume:<5f}\n")
        else : summary_file.write(f"-> Summed volume has not been calculated yet\n")
        summary_file.write(f"elementary volume                          {self.d_V:<5f}\n")
        summary_file.write(f"total_electrons (from siesta)              {self.number_of_electrons}\n")
        summary_file.write(f"normalized total electronic charge         {self.total_normalized_electronic_charge}\n")
        summary_file.write(f"total ionic charge                         {self.summed_ionic_charge}\n")
        summary_file.write(f"Normalization factor for electron density  {self.charge_normalization_factor:<5f}\n")
        summary_file.write(f"Total unnormalized charge from cube file   {self.total_unnormalized_charge:<5f}\n")
        summary_file.write(f"Total Charge                               {self.integrated_charge}\n")
        summary_file.write(f"dipole                                     {self.dipole*unit_factor_Debye} Debye\n")
        summary_file.write(f"dipole                                     {self.dipole} in q.r numofelectrons.angs\n")
        summary_file.write(f"electronic dipole                          {self.electronic_dipole*unit_factor_Debye} Debye\n")
        summary_file.write(f"electronic dipole                          {self.electronic_dipole} in q.r numofelectrons.angs\n")
        summary_file.write(f"ionic dipole                               {self.ionic_dipole*unit_factor_Debye} Debye\n")
        summary_file.write(f"ionic dipole                               {self.ionic_dipole} in q.r numofelectrons.angs\n")
        for i,(k,v)in enumerate(conversions.items()):
            summary_file.write(f"\n")
            summary_file.write(f"Q                              in {k}\n")
            summary_file.write(f"{self.Q*v}\n")
            summary_file.write(f"\n")
            summary_file.write(f"Q (in the non-traceless from)  in {k}\n")
            summary_file.write(f"{self.Q_non_traceless*v}\n")
            summary_file.write(f"\n")
            for col in range(0,3):
                for row in range(0,3):
                    summary_file.write(f"Q{get_direction(row)}{get_direction(col)}                                        {self.Q[row,col]*v:<5f} {k}\n")
            for col in range(0,3):
                for row in range(0,3):
                    summary_file.write(f"Q{get_direction(row)}{get_direction(col)}  (in the non-traceless from)           {self.Q_non_traceless[row,col]*v:<5f} {k}\n")
        summary_file.close()

        # Writing out the new data in a new file
        with open(f"{self.out_file_name}.electronicandionic.cube", "w+") as file:
            file.write("Produced by Eranjan\n")
            file.write(f"For total charge for system: {self.out_file_name}\n")
            for line in self.first_lines_of_file:
                file.write(line)
            for x in range(self.a_number_of_voxels):
                for y in range(self.b_number_of_voxels):
                    for z in range(self.c_number_of_voxels):            
                        file.write(f"{self.rho[x][y][z]:1.5E} ")
                        if (z % 6 == 5):
                         file.write("\n")
                    file.write("\n")

    def get_r_vec(self,x,y,z):
        """"""
        r_vec = x*self.a_voxel_vec + y*self.b_voxel_vec + z*self.c_voxel_vec + self.origin
        return r_vec

    def calculate_volume(self):
        """
        This method sums the volume from using the volume of a single voxel and iterating through the volxels.
        """
        self.summed_volume=0
        prog = progress_bar(self.a_number_of_voxels*self.b_number_of_voxels*self.c_number_of_voxels, descriptor="Calculating Summed volume")
        for ia in range(self.a_number_of_voxels):
            for ib in range(self.b_number_of_voxels):
                for ic in range(self.c_number_of_voxels):
                    prog.get_progress((ia)*(self.b_number_of_voxels*self.c_number_of_voxels)+(ib)*(self.c_number_of_voxels)+ic)
                    self.summed_volume+=self.d_V 
        return self.summed_volume
                    
    def calculate_center_of_electronic_charge(self):
        """"""
        self.electronic_dipole = 0

        prog = progress_bar(self.a_number_of_voxels*self.b_number_of_voxels*self.c_number_of_voxels, descriptor="Calculating Centre of Charge")
        for ia in range(self.a_number_of_voxels):
            for ib in range(self.b_number_of_voxels):
                for ic in range(self.c_number_of_voxels):
                    prog.get_progress((ia)*(self.b_number_of_voxels*self.c_number_of_voxels)+(ib)*(self.c_number_of_voxels)+ic)
                    r_vec = self.get_r_vec(ia, ib, ic)
                    self.electronic_dipole+=self.rho[ia, ib, ic]*r_vec
        self.center_of_charge_electronic = self.electronic_dipole/self.total_normalized_electronic_charge
        return self.center_of_charge_electronic

    def calculate_dipole_moments(self):
        """"""
        self.dipole = 0
        self.total_normalized_electronic_charge = 0

        prog = progress_bar(self.a_number_of_voxels*self.b_number_of_voxels*self.c_number_of_voxels, descriptor="Calculating dipole moments")
        for ia in range(self.a_number_of_voxels):
            for ib in range(self.b_number_of_voxels):
                for ic in range(self.c_number_of_voxels):
                    prog.get_progress((ia)*(self.b_number_of_voxels*self.c_number_of_voxels)+(ib)*(self.c_number_of_voxels)+ic)
                    r_vec = self.get_r_vec(ia, ib, ic)
                    self.dipole+=self.rho[ia, ib, ic]*r_vec

        # Now adding the ionic part
        for i_atom, atom in enumerate(self.atoms_info):
            coordinates = atom[3]
            charge = atom[2]   # This is from the valance electron item in the list
            r_vec = np.array([coordinates[0], coordinates[1], coordinates[2]])
            self.dipole += charge*r_vec
        return self.dipole

    def calculate_quadrupole_moments(self):
        """"""
        self.Q = np.zeros((3,3))
        self.Q_non_traceless = np.zeros((3,3))
        self.dipole = np.zeros((1,3))
        self.calculate_electronic_part()
        self.calculate_ionic_part()
        

    def calculate_electronic_part(self):
        # This part calculates only the electronic part from the cube file
        self.electronic_dipole = np.zeros((1,3))
        prog = progress_bar(self.a_number_of_voxels*self.b_number_of_voxels*self.c_number_of_voxels, descriptor="Analyzing for moments")
        for ia in range(self.a_number_of_voxels):
            for ib in range(self.b_number_of_voxels):
                for ic in range(self.c_number_of_voxels):
                    prog.get_progress((ia)*(self.b_number_of_voxels*self.c_number_of_voxels)+(ib)*(self.c_number_of_voxels)+ic)
                    """Methana podi indeces prashnayak thiyeanwa. mokenda iterate karanna one i+1 da nattam i da kiyala"""
                    r_vec = self.get_r_vec(ia, ib, ic)
                    r2 = np.dot(r_vec, r_vec)
                    r = np.sqrt(r2)
                    self.dipole+=self.rho[ia, ib, ic]*r_vec*self.d_V
                    self.electronic_dipole+=self.rho[ia, ib, ic]*r_vec*self.d_V
                    for col in range(0,3):
                        for row in range(0,3):
                            if col == row: f=1
                            else:f=0
                            self.Q[row,col]               += self.rho[ia, ib, ic]*(3*(r_vec[row])*(r_vec[col]) - r2*f)*self.d_V
                            self.Q_non_traceless[row,col] += self.rho[ia, ib, ic]*((r_vec[row])*(r_vec[col]))*self.d_V

    def calculate_ionic_part(self):
        # Now adding the ionic part
        self.summed_ionic_charge = 0
        self.ionic_dipole = np.zeros((1,3))
        for i_atom, atom in enumerate(self.atoms_info):
            coordinates = atom[3]
            charge = atom[2]   # This is from the valance electron item in the list
            r_vec = np.array([coordinates[0], coordinates[1], coordinates[2]])
            r2 = np.dot(r_vec,r_vec)
            self.dipole += charge*r_vec
            self.ionic_dipole += charge*r_vec
            self.summed_ionic_charge += charge
            for col in range(0,3):
                for row in range(0,3):
                    if col == row: f=1
                    else:f=0
                    self.Q[row,col]               += charge*(3*(r_vec[row])*(r_vec[col]) - r2*f)
                    self.Q_non_traceless[row,col] += charge*((r_vec[row])*(r_vec[col]))

    def get_all_moments(self, file_name = None, out_put_file_name = "", cell = [[0., 0., 0.,],[0., 0., 0.,],[0., 0., 0.,]]):
        """
        This is the main method to calculate the dipole/quadrupole moments
        """
        import numpy as np

        print("\nCalculating the dipoles and quadrupoles...")
        self.read_in_rho_cube_file(file_name=file_name)
        self.normalize_electronic_charge_density()
        self.get_volumes_from_cube_header()
        self.calculate_volume()
        self.calculate_quadrupole_moments()
        self.integrated_charge = self.summed_ionic_charge + self.total_normalized_electronic_charge
        self.write_to_quadrupole_outputfiles(out_put_file_name)

        return_dict = {}
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

        self.Q = np.zeros((3,3))
        self.Q_non_traceless = np.zeros((3,3))

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
            if "electronic dipole" in line and "Debye" in line and "unadjusted" not in line and "binned" not in line:
                line = line.strip("[")
                line = line.strip("]")
                # print(line)
                parsed = line.split()
                if "[" in parsed : parsed.remove("[")
                if "]" in parsed : parsed.remove("]")
                # print(parsed)
                self.electronic_dipole = [ float(parsed[-4].strip("[")), float(parsed[-3]), float(parsed[-2].strip("]")) ]
            if "ionic dipole" in line and "Debye" in line and "unadjusted" not in line and "binned" not in line:
                line = line.strip("[")
                line = line.strip("]")
                # print(line)
                parsed = line.split()
                if "[" in parsed : parsed.remove("[")
                if "]" in parsed : parsed.remove("]")
                # print(parsed)
                self.ionic_dipole = [ float(parsed[-4].strip("[")), float(parsed[-3]), float(parsed[-2].strip("]")) ]
            if "Qxx" in line and "Debye.Angs" in line and "(in the non-traceless from)" not in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Q_non_traceless[0,0] = float(parsed[1])
            if "Qyy" in line and "Debye.Angs" in line and "(in the non-traceless from)" not in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Q_non_traceless[1,1] = float(parsed[1])
            if "Qzz" in line and "Debye.Angs" in line and "(in the non-traceless from)" not in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Q_non_traceless[2,2] = float(parsed[1])
            if "Qxy" in line and "Debye.Angs" in line and "(in the non-traceless from)" not in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Q_non_traceless[0,1] = float(parsed[1])
            if "Qxz" in line and "Debye.Angs" in line and "(in the non-traceless from)" not in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Q_non_traceless[0,2] = float(parsed[1])
            if "Qyx" in line and "Debye.Angs" in line and "(in the non-traceless from)" not in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Q_non_traceless[1,0] = float(parsed[1])
            if "Qyz" in line and "Debye.Angs" in line and "(in the non-traceless from)" not in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Q_non_traceless[1,2] = float(parsed[1])
            if "Qzx" in line and "Debye.Angs" in line and "(in the non-traceless from)" not in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Q_non_traceless[2,0] = float(parsed[1])
            if "Qzy" in line and "Debye.Angs" in line and "(in the non-traceless from)" not in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Q_non_traceless[2,1] = float(parsed[1])


            if "Qxx" in line and "Debye.Angs" in line and "(in the non-traceless from)" in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Q_non_traceless[0,0] = float(parsed[5])
            if "Qyy" in line and "Debye.Angs" in line and "(in the non-traceless from)" in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Q_non_traceless[1,1] = float(parsed[5])
            if "Qzz" in line and "Debye.Angs" in line and "(in the non-traceless from)" in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Q_non_traceless[2,2] = float(parsed[5])
            if "Qxy" in line and "Debye.Angs" in line and "(in the non-traceless from)" in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Q_non_traceless[0,1] = float(parsed[5])
            if "Qxz" in line and "Debye.Angs" in line and "(in the non-traceless from)" in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Q_non_traceless[0,2] = float(parsed[5])
            if "Qyx" in line and "Debye.Angs" in line and "(in the non-traceless from)" in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Q_non_traceless[1,0] = float(parsed[5])
            if "Qyz" in line and "Debye.Angs" in line and "(in the non-traceless from)" in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Q_non_traceless[1,2] = float(parsed[5])
            if "Qzx" in line and "Debye.Angs" in line and "(in the non-traceless from)" in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Q_non_traceless[2,0] = float(parsed[5])
            if "Qzy" in line and "Debye.Angs" in line and "(in the non-traceless from)" in line and "numofelectrons.angs^2" not in line:
                parsed = line.split()
                self.Q_non_traceless[2,1] = float(parsed[5])
            if "Total unnormalized charge from cube file" in line:
                parsed = line.split()
                self.total_unnormalized_charge = float(parsed[-1])
            if "Full volume (summed)" in line:
                parsed = line.split()
                self.summed_volume = float(parsed[-1])

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
