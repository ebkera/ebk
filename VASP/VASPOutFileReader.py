"""
This file reads the the file system_label.out and extracts/calculates data/values from it
"""

from ebk.BandPlotter import BandPlotter 
import numpy as np

class VASPReadOut():
    def __init__(self, out_folder, *args, **kwargs):
        self.out_folder = out_folder
        self.bands = []
        self.kpoints = []
        trigger = False
        lattice_vectors_trigger = False
        self.k_dist=[0]
        self.highest_valance = [0, -500]    # Will contain [[kpath_index],E]
        self.lowest_conduction = [0, 500]  # Will contain [[kpath_index],E]
        self.valance_band_label_index = 0
        self.conduction_band_label_index = 0
        self.hsp = []
        self.hss = []
        self.Eg = "notset"
        self.draw_band_edge_lines = False
        self.spin_orbit = None
        self.Ef = 0
        self.E_tot = 0
        self.E_tot_with_entropy = 0
        self.consider_actual_k_distance = False  # important when we want actual K distance like when calculating electron and hole masses
        self.kpoints_read_from = kwargs.get("kpoints_read_from", 1)  # This is when you want to ignore some of the first k points for exxample in a HSE calculation grid k points not requried for band structure calculations
        self.reset_k_dist = kwargs.get("reset_k_dist", [])   # This is when there is a discontinuity in the K path like U|K
        self.lattice_vectors = []
        self.reciprocal_lattice_vectors = []
        
        with open(f"{self.out_folder}/OUTCAR", "r+") as OUTCAR:
            for line in OUTCAR:
                # Getting paramters:
                if "LSORBIT" in line:
                    k = line.split()
                    if k[2] == "T":
                        self.spin_orbit = True
                    elif k[2] == "F":
                        self.spin_orbit = False
                    else:
                        print(f"Spin-Orbit settings not detected")

                # Reading in the lattice vectors
                # a1,a2,a3 are the real vectors in 2pi/a while A1,A2,A3 are the real cell lengths in Angstroms
                if 'A1 = ' in line:
                    a = line.split()
                    self.a1= [float(a[3].replace(',',"")), float(a[4].replace(',',"")), float(a[5].replace(')',""))]
                    self.A1 = np.sqrt(self.a1[0]**2+self.a1[1]**2+self.a1[2]**2)
                if 'A2 = ' in line:
                    a = line.split()
                    self.a2= [float(a[3].replace(',',"")), float(a[4].replace(',',"")), float(a[5].replace(')',""))]
                    self.A2 = np.sqrt(self.a2[0]**2+self.a2[1]**2+self.a2[2]**2)
                if 'A3 = ' in line:
                    a = line.split()
                    self.a3= [float(a[3].replace(',',"")), float(a[4].replace(',',"")), float(a[5].replace(')',""))]
                    self.A3 = np.sqrt(self.a3[0]**2+self.a3[1]**2+self.a3[2]**2)

                # This part will read in the lines from the initial k point listing where we also have the plane wave numbers given
                # Sometimes VASP output is not properly syntaxed and there are not spaces for so will not be able to resolve values
                # So this has been abandoned and we read K points from the second section with energies
                # if 'k-point' in line and 'plane waves' in line:
                #     k = line.split()
                #     k2 = [float(k[3]),float(k[4]),float(k[5])]
                #     self.kpoints.append(k2)
                #     # print(k2)

                # print(line)  # Left here for debugging 
                if 'k-point' in line and len(line.split()) == 6 and "Found" not in line  and "KPOINTS" not in line and "generate" not in line:
                    k = line.split()
                    current_kpoint = int(k[1])
                    if current_kpoint< self.kpoints_read_from : continue
                    k2 = [float(k[3]),float(k[4]),float(k[5])]
                    self.kpoints.append(k2)
                    # print(k2)

                if 'E-fermi' in line:
                    k = line.split()
                    E_fermi = float(k[2])
                    self.Ef = E_fermi
                    # E_fermi = 0
                    # print(E_fermi)
                    
                if 'TOTEN  =' in line:
                    k = line.split()
                    self.E_tot_with_entropy = float(k[4])
                    
                if 'energy  without entropy=' in line:
                    k = line.split()
                    self.E_tot = float(k[3])

                if f" k-point     {self.kpoints_read_from} " in line:
                    # print(line)
                    trigger = True
                    soft_trigger = True
                    
                if "--------------------------------------------------------------------------------------------------------" in line and trigger == True:
                    trigger = False
                    # print(line)

                if trigger == True:
                    if "k-point" in line:
                        band_count = 0

                if trigger == True and len(line.split()) == 3:
                    f=line.split()
                    # print(f)
                    band_label = int(f[0])-1
                    # E = float(f[1])-E_fermi
                    E = float(f[1])
                    # E = float(f[1])  # if we want to not shift by Ef
                    Occ = float(f[2])
                    try:
                        self.bands[band_label].append(E)
                    except:
                        self.bands.append([E])
                    
                    if self.spin_orbit:
                        Occ_threshold = 0.9
                    elif self.spin_orbit == False:
                        Occ_threshold = 1.9
                    else: 
                        # For cases where SO is not detected
                        Occ_threshold = 1.9

                    if Occ > Occ_threshold:  # This is because there might be partial occupancies
                        if E >= self.highest_valance[1]: 
                            self.highest_valance[1] =  E
                            self.highest_valance[0] = current_kpoint - 1
                            self.valance_band_label_index = band_label
                    else: # This is because there might be partial occupancies
                        if E<= self.lowest_conduction[1]:
                            self.lowest_conduction[1] = E
                            self.lowest_conduction[0] = current_kpoint - 1
                            self.conduction_band_label_index = band_label
                            # print(E)

                if "direct lattice vectors" in line and "reciprocal lattice vectors" in line: 
                    lattice_vectors_trigger = True
                    self.lattice_vectors = []
                if lattice_vectors_trigger == True and "direct lattice vectors" not in line:
                    parsed = line.split()
                    self.lattice_vectors.append([float(parsed[0]), float(parsed[1]), float(parsed[2])])
                    self.reciprocal_lattice_vectors.append([float(parsed[3]), float(parsed[4]), float(parsed[5])])
                    if len(self.lattice_vectors) == 3: lattice_vectors_trigger = False
                
                # This part was started to write the optical data from an optical run into seperate files, but was abandoned and has its own function now.
                # if "frequency dependent IMAGINARY DIELECTRIC FUNCTION" in line:
                #     REAL_trigger = False
                #     IMAG_trigger = True
                # elif "frequency dependent      REAL DIELECTRIC FUNCTION" in line:
                #     REAL_trigger = True
                #     IMAG_trigger = False
                # if REAL_trigger:
                #     if line != "": er_data.append(line)
                #     else:REAL_trigger = False
                # if IMAG_trigger:
                #     if line != "": ei_data.append(line)
                #     else: IMAG_trigger = False
              
        self.Eg = self.lowest_conduction[1] - self.highest_valance[1]

        def process_symbols(hss):
            if hss.lower() == "gamma" or hss.lower() == "g": return "$\Gamma$"
            else: return hss.upper()
        
        def get_k_dist():
           # Here we adjust for the actual lengths in reciprocal space
           # While the correct thing would be to leave everything in terms of reciprocal space this comes in handy when trying to compare band structure between software
            if self.consider_actual_k_distance:
                x_multiplier = 2*np.pi**(10)/self.A1
                y_multiplier = 2*np.pi**(10)/self.A2
                z_multiplier = 2*np.pi**(10)/self.A3
            else:
                x_multiplier = 1
                y_multiplier = 1
                z_multiplier = 1
            dist_v=[0, 0, 0]
            k_start = (len(self.hss)-2)*self.k_point_density
            for x in range(k_start, k_start+self.k_point_density):
                if x == 0: continue
                for i in range(3):
                    dist_v[i] = self.kpoints[x][i] - self.kpoints[x-1][i]
                if x == k_start: self.k_dist.append(self.k_dist[-1]) # here discontinuities are taken into account
                else: self.k_dist.append(((dist_v[0]*x_multiplier)**2+(dist_v[1]*y_multiplier)**2+(dist_v[2]*z_multiplier)**2) + self.k_dist[x-1])
           
        try:
            with open(f"{self.out_folder}/KPOINTS", "r+") as KPOINTS:
                hss_counter = 0
                for line in KPOINTS:
                    line_s = line.split()
                    if len(line_s) == 1:
                        try: self.k_point_density = int(line_s[0])
                        except: pass
                    try: float(line_s[0])  # This will check for actual K points
                    except: continue
                    
                    if len(line_s) == 4:
                        hss_counter +=1
                        hss = line_s[3]
                        if hss_counter == 1: 
                            self.hss.append(process_symbols(hss))
                            self.hsp.append(0) 
                        elif hss_counter % 2 == 0: # this is when a single direction is closed between two points happen
                            self.hss.append(process_symbols(hss))
                            get_k_dist()
                            self.hsp.append(self.k_dist[-1])
                        elif hss_counter % 2 == 1: # this is beginign of a new set of two points
                            if process_symbols(hss) != self.hss[-1] : # here is seems there is a discontinuity
                                self.hss[-1] = f"{self.hss[-1]} | {process_symbols(hss)}"

        except FileNotFoundError as e:
            print(f"Could not find KPOINTS file for : {self.out_folder}/KPOINTS")
            print(f"{e}")

    def get_dos(self, DOS_DIR):
        self.dos_MIN_E = 0
        self.dos_MAX_E = 0
        self.dos_E = []
        self.dos_D = []  # These are the DOS values
        self.dos_ID = []  # These are the integrated DOS values
        # print("Calculating Density of States")
        with open(f"{DOS_DIR}/DOSCAR", "r+") as DOSCAR:
            for i,line in enumerate(DOSCAR):
                if i == 5:
                    data = line.split()
                    self.dos_MIN_E = float(data[0])
                    self.dos_MAX_E = float(data[1])
                    self.dos_N_points = int(data[2])
                    self.Ef_dos = float(data[3]) # Ef is updated here
                if i>5 and len(line.split()) == 3:
                    data = line.split()
                    self.dos_E.append(float(data[0]) - self.Ef_dos)
                    self.dos_D.append(float(data[1]))
                    self.dos_ID.append(float(data[2]))

    def plot_dos(self):
        import matplotlib.pyplot as plt
        plt.plot(self.dos_E, self.dos_D)
        plt.show()

    def get_final_volume(self):
        import numpy as np
        return np.dot(self.lattice_vectors[2], np.cross(self.lattice_vectors[0],self.lattice_vectors[1]))

    def convert_to_cif(self, label = "structure"):
        from ase.io import read, write
        structure =  read(f"{self.out_folder}/CONTCAR")
        write(f"{self.out_folder}/{label}.cif", structure)

    def get_band_gap(self):
        return self.Eg

    def get_fermi_energy(self):
        return self.Ef

    def get_highest_valance_band_kpoint(self):
        return self.highest_valance[0]

    def get_lowest_conduction_band_kpoint(self):
        return self.lowest_conduction[0]

    def get_highest_valance_band_energy(self):
        return self.highest_valance[1]

    def get_lowest_conduction_band_energy(self):
        return self.lowest_conduction[1]

    def get_total_energy(self):
        return self.E_tot

    def get_total_energy_without_entropy(self):
        return self.E_tot_with_entropy

    def get_cell_volume(self):
        return np.dot(self.a1, np.cross(self.a2, self.a3) )

    def get_band_structure(self, file_name = None):
        """
        Sets initial settings and returns a BandPlotter type object for further manipulation but ready to plot
        """
        if file_name == None:
            file_name = self.out_folder
        x = BandPlotter()
        if hasattr(self, "dos_E"):
            x.add_to_plot(self.highest_valance[1], self.k_dist, self.bands, self.hsp, self.hss, self.dos_E, self.dos_D)
            # print("methani weda karanne")
            # print(self.dos_D)
        else:
            if self.Eg >= 0: 
                x.add_to_plot(self.Ef, self.k_dist, self.bands, self.hsp, self.hss, label=f"$E_{{g}}={self.get_band_gap(): 3.3f} eV$")
                arrow_data = [0,self.get_highest_valance_band_energy()-self.Ef, self.k_dist[self.get_lowest_conduction_band_kpoint()], self.get_lowest_conduction_band_energy()-self.Ef]
            else:
                x.add_to_plot(self.Ef, self.k_dist, self.bands, self.hsp, self.hss, label=f"$E_{{g}}={0: 3.3f} eV)$")
                arrow_data = [0,0, 0, 0]
            x.arrow_data.append(arrow_data)
        x.file_name = f"{file_name}"
        x.set_y_range=True
        x.same_band_colour = True
        x.hlines = True
        x.vlines = True
        x.ylim_high = 2.5
        x.ylim_low = -2.5
        x.saveas_extension = "png"
        if self.draw_band_edge_lines:
            x.extra_horizontal_lines = [[self.lowest_conduction[1], f"Conduction band edge"]]
            x.extra_horizontal_lines.append([self.highest_valance[1], f"Valance band edge"])
        return x

    def get_curvature_at_k_point(self, k_point: list[float,float,float]=False, k_index:int=False, direction:str=False, sigma:int=15, high_symmetry_point:bool=True, show_plot:bool=False, show_holes:bool=False, show_elec:bool=False)-> float:
        """
        This code is still being written.
        Goals: to get the band curvature around the maxima and minima of the CBM and VBM
        k_point: 
        direction: A string of length 2
        high_symmetry_point: set to true to accomodate for VASP k points being repeated in the k-path at highsymmetry points
        """
        print(f"Calculating curvature of the bands")
        import numpy as np
        import scipy.constants

        if k_point:
            pass

        if k_index:
            print(f"Calculating around k point: {self.kpoints[k_index]}")
            print(f"Please check that valance and conduction band indeces are: {self.valance_band_label_index+1, self.conduction_band_label_index+1}")
            right_side_conduction_E = [self.bands[self.conduction_band_label_index][x] for x in range(k_index, k_index+sigma)]
            right_side_conduction_E_ori = right_side_conduction_E.copy()
            right_side_kpoint_dist = [self.k_dist[x] for x in range(k_index, k_index+sigma)]
            right_side_kpoint_dist_delta = right_side_kpoint_dist[-1] - right_side_kpoint_dist[-2]
            # prepending
            for x in range(1,len(right_side_conduction_E)):
                right_side_conduction_E.insert(0,right_side_conduction_E_ori[x])
                right_side_kpoint_dist.insert(0,right_side_kpoint_dist[0]-right_side_kpoint_dist_delta)
            if high_symmetry_point:
                left_side_conduction_E = [self.bands[self.conduction_band_label_index][x] for x in range(k_index-sigma, k_index)]
                left_side_conduction_E_ori = left_side_conduction_E.copy()
                left_side_kpoint_dist = [self.k_dist[x] for x in range(k_index-sigma, k_index)]
            else:
                left_side_conduction_E = [self.bands[self.conduction_band_label_index][x] for x in range(k_index-sigma+1, k_index+1)]
                left_side_conduction_E_ori = left_side_conduction_E.copy()
                left_side_kpoint_dist = [self.k_dist[x] for x in range(k_index-sigma+1, k_index+1)]
            left_side_kpoint_dist_delta = left_side_kpoint_dist[-1] - left_side_kpoint_dist[-2]
            #appending
            for x in range(len(left_side_conduction_E_ori)-2, -1, -1):
                left_side_conduction_E.append(left_side_conduction_E_ori[x])
                left_side_kpoint_dist.append(left_side_kpoint_dist[-1]+left_side_kpoint_dist_delta)

            right_side_valance_E = [self.bands[self.valance_band_label_index][x] for x in range(k_index, k_index+sigma)]
            right_side_valance_E_ori = right_side_valance_E.copy()
            # prepending
            for x in range(1,len(right_side_valance_E)):
                right_side_valance_E.insert(0,right_side_valance_E_ori[x])
            if high_symmetry_point:
                left_side_valance_E = [self.bands[self.valance_band_label_index][x] for x in range(k_index-sigma, k_index)]
                left_side_valance_E_ori = left_side_valance_E.copy()
            else:
                left_side_valance_E = [self.bands[self.valance_band_label_index][x] for x in range(k_index-sigma+1, k_index+1)]
                left_side_valance_E_ori = left_side_valance_E.copy()
            #appending
            for x in range(len(left_side_valance_E_ori)-2, -1, -1):
                left_side_valance_E.append(left_side_valance_E_ori[x])

            fit_para_left_conduction = np.polyfit(left_side_kpoint_dist, left_side_conduction_E, 2)
            fit_para_right_conduction = np.polyfit(right_side_kpoint_dist, right_side_conduction_E, 2)
            fit_para_left_valance = np.polyfit(left_side_kpoint_dist, left_side_valance_E, 2)
            fit_para_right_valance = np.polyfit(right_side_kpoint_dist, right_side_valance_E, 2)

            electron_effective_mass_left = scipy.constants.hbar**2/(2*fit_para_left_conduction[0]*scipy.constants.elementary_charge*scipy.constants.m_e)  # The two on the denominator is from the 2 that comes from differentiating a 2 order polynomial
            electron_effective_mass_right = scipy.constants.hbar**2/(2*fit_para_right_conduction[0]*scipy.constants.elementary_charge*scipy.constants.m_e)  # The two on the denominator is from the 2 that comes from differentiating a 2 order polynomial
            hole_effective_mass_left = scipy.constants.hbar**2/(2*fit_para_left_valance[0]*scipy.constants.elementary_charge*scipy.constants.m_e)  # The two on the denominator is from the 2 that comes from differentiating a 2 order polynomial
            hole_effective_mass_right = scipy.constants.hbar**2/(2*fit_para_right_valance[0]*scipy.constants.elementary_charge*scipy.constants.m_e)  # The two on the denominator is from the 2 that comes from differentiating a 2 order polynomial

            # electron_effective_mass_left = scipy.constants.hbar**2/(2*fit_para_left_conduction[0])  # The two on the denominator is from the 2 that comes from differentiating a 2 order polynomial
            # electron_effective_mass_right = scipy.constants.hbar**2/(2*fit_para_right_conduction[0])  # The two on the denominator is from the 2 that comes from differentiating a 2 order polynomial
            # hole_effective_mass_left = scipy.constants.hbar**2/(2*fit_para_left_valance[0])  # The two on the denominator is from the 2 that comes from differentiating a 2 order polynomial
            # hole_effective_mass_right = scipy.constants.hbar**2/(2*fit_para_right_valance[0])  # The two on the denominator is from the 2 that comes from differentiating a 2 order polynomial

            return_dict = {}
            return_dict.update({"fit_para_conduction_l":fit_para_left_conduction})
            return_dict.update({"fit_para_conduction_r":fit_para_right_conduction})
            return_dict.update({"fit_para_valance_l":fit_para_left_valance})
            return_dict.update({"fit_para_valance_r":fit_para_left_valance})
            return_dict.update({"me_l":electron_effective_mass_left})
            return_dict.update({"me_r":electron_effective_mass_right})
            return_dict.update({"mh_l":hole_effective_mass_left})
            return_dict.update({"mh_r":hole_effective_mass_right})

            p_era_LDA_ES = np.poly1d(fit_para_left_conduction)
            fit_curve_left_conduction = p_era_LDA_ES(left_side_kpoint_dist)
            p_right_conduction = np.poly1d(fit_para_right_conduction)
            fit_curve_right_conduction = p_right_conduction(right_side_kpoint_dist)
            
            p_left_valance = np.poly1d(fit_para_left_valance)
            fit_curve_left_valance = p_left_valance(left_side_kpoint_dist)
            p_right_valance = np.poly1d(fit_para_right_valance)
            fit_curve_right_valance = p_right_valance(right_side_kpoint_dist)
            import matplotlib.pyplot as plt 
            if show_elec:
                plt.plot(left_side_kpoint_dist, fit_curve_left_conduction, label=f"Left Conduction fit")
                plt.plot(left_side_kpoint_dist, left_side_conduction_E, label=f"Left Conduction")
                plt.plot(right_side_kpoint_dist, fit_curve_right_conduction, label=f"Right Conduction fit")
                plt.plot(right_side_kpoint_dist, right_side_conduction_E, label=f"Right Conduction")
            if show_holes:
                plt.plot(left_side_kpoint_dist, fit_curve_left_valance, label=f"Left Valance fit")
                plt.plot(left_side_kpoint_dist, left_side_valance_E, label=f"Left Valance")
                plt.plot(right_side_kpoint_dist, fit_curve_right_valance, label=f"Right Valance fit")
                plt.plot(right_side_kpoint_dist, right_side_valance_E, label=f"Right Valance")
            plt.legend()
            if show_plot: plt.show()

            return return_dict

        if direction:
            pass

    def get_optical(self, OPTICS_DIR):
        from ebk.BandPlotter import BandPlotter
        import matplotlib.pyplot as plt
        # E_1 is the real part of the dielectric
        # E_2 is the imaginary part of the dielectric
        # E is energy
        print("Reading in Optics data...")
        E_1_x = []  
        E_2_x = []
        E_1_y = []
        E_2_y = []
        E_1_z = []
        E_2_z = []
        E_1_xy = []
        E_2_xy = []
        E_1_yz = []
        E_2_yz = []
        E_1_zx = []
        E_2_zx = []
        E = []

        with open(f"{OPTICS_DIR}/OUTCAR", "r+") as OUTCAR:
            for i,line in enumerate(OUTCAR):
                REAL_trigger = False
                IMAG_trigger = False
                new_trigger  = False
                E_1_x_temp = []
                E_2_x_temp = []
                E_1_y_temp = []
                E_2_y_temp = []
                E_1_z_temp = []
                E_2_z_temp = []
                E_1_xy_temp = []
                E_2_xy_temp = []
                E_1_yz_temp = []
                E_2_yz_temp = []
                E_1_zx_temp = []
                E_2_zx_temp = []
                E_temp = []

                for line in OUTCAR:
                    if "frequency dependent IMAGINARY DIELECTRIC FUNCTION" in line:
                        REAL_trigger = False
                        IMAG_trigger = True
                        new_trigger  = True
                        continue
                    elif "frequency dependent      REAL DIELECTRIC FUNCTION" in line:
                        REAL_trigger = True
                        IMAG_trigger = False
                        new_trigger  = True
                        continue           
                    if REAL_trigger or IMAG_trigger:
                        try:
                            split_line = line.split()
                            temp_line_vals = []
                            if len(split_line) == 0:
                                if REAL_trigger:
                                    E_1_x.append(E_1_x_temp)
                                    E_1_y.append(E_1_y_temp)
                                    E_1_z.append(E_1_z_temp)
                                    E_1_xy.append(E_1_xy_temp)
                                    E_1_yz.append(E_1_yz_temp)
                                    E_1_zx.append(E_1_zx_temp)
                                if IMAG_trigger:
                                    E_2_x.append(E_2_x_temp)
                                    E_2_y.append(E_2_y_temp)
                                    E_2_z.append(E_2_z_temp)
                                    E_2_xy.append(E_2_xy_temp)
                                    E_2_yz.append(E_2_yz_temp)
                                    E_2_zx.append(E_2_zx_temp)
                                E.append(E_temp)
                                REAL_trigger = False
                                IMAG_trigger = False
                                E_1_x_temp = []
                                E_2_x_temp = []
                                E_1_y_temp = []
                                E_2_y_temp = []
                                E_1_z_temp = []
                                E_2_z_temp = []
                                E_1_xy_temp = []
                                E_2_xy_temp = []
                                E_1_yz_temp = []
                                E_2_yz_temp = []
                                E_1_zx_temp = []
                                E_2_zx_temp = []
                                E_temp = []

                            if len(split_line) != 7: continue
                            for x in split_line:
                                # print("working")
                                temp_line_vals.append(float(x))
                            # print(temp_line_vals)
                        except: continue

                        if REAL_trigger:                    
                            E_1_x_temp.append(temp_line_vals[1])
                            E_1_y_temp.append(temp_line_vals[2])
                            E_1_z_temp.append(temp_line_vals[3])
                            E_1_xy_temp.append(temp_line_vals[4])
                            E_1_yz_temp.append(temp_line_vals[5])
                            E_1_zx_temp.append(temp_line_vals[6])
                        if IMAG_trigger:                    
                            E_2_x_temp.append(temp_line_vals[1])
                            E_2_y_temp.append(temp_line_vals[2])
                            E_2_z_temp.append(temp_line_vals[3])
                            E_2_xy_temp.append(temp_line_vals[4])
                            E_2_yz_temp.append(temp_line_vals[5])
                            E_2_zx_temp.append(temp_line_vals[6])
                        E_temp.append(temp_line_vals[0])

        plt.figure()
        plt.plot(E[0], E_2_x[0], label=f"$\epsilon_2 XX$ (d-d)")
        plt.plot(E[0], E_2_y[0], label=f"$\epsilon_2 YY$ (d-d)")
        plt.plot(E[0], E_2_z[0], label=f"$\epsilon_2 ZZ$ (d-d)")
        # plt.plot(E[-1], E_2_x[1], label=f"$\epsilon_2 XX$ (c-c)")
        # plt.plot(E[-1], E_2_y[1], label=f"$\epsilon_2 YY$ (c-c)")
        # plt.plot(E[-1], E_2_z[1], label=f"$\epsilon_2 ZZ$ (c-c)")
        plt.legend()
        # plt.title(f"$\epsilon_2$ vs Energy")
        plt.xlabel("Energy (eV)")
        plt.ylabel("$\epsilon_2$")
        plt.xlim([0, 5])
        plt.savefig(f"{OPTICS_DIR}/E2vsE.png")
        plt.show()

        plt.figure()
        plt.plot(E[0], E_1_x[0], label=f"$\epsilon_1 XX$ (d-d)")
        plt.plot(E[0], E_1_y[0], label=f"$\epsilon_1 YY$ (d-d)")
        plt.plot(E[0], E_1_z[0], label=f"$\epsilon_1 ZZ$ (d-d)")
        # plt.plot(E[-1], E_1_x[1], label=f"$\epsilon_1 XX$ (c-c)")
        # plt.plot(E[-1], E_1_y[1], label=f"$\epsilon_1 YY$ (c-c)")
        # plt.plot(E[-1], E_1_z[1], label=f"$\epsilon_1 ZZ$ (c-c)")
        plt.legend()
        # plt.title(f"$\epsilon_1$ vs Energy")
        plt.xlabel("Energy (eV)")
        plt.ylabel("$\epsilon_1$")
        plt.xlim([0, 5])
        plt.savefig(f"{OPTICS_DIR}/E1vsE.png")
        # plt.show()

        with open(f"{OPTICS_DIR}/e_r.dat", "w+") as er:
            for i,v in enumerate(E[0]):
                if i !=  0:er.write("\n")
                er.write(f"{E[0][i]:>12.6f}{E_1_x[0][i]:>12.6f}{E_1_y[0][i]:>12.6f}{E_1_z[0][i]:>12.6f}{E_1_xy[0][i]:>12.6f}{E_1_yz[0][i]:>12.6f}{E_1_zx[0][i]:>12.6f}")

        with open(f"{OPTICS_DIR}/e_i.dat", "w+") as ei:
            for i,v in enumerate(E[0]):
                if i != 0: ei.write("\n")
                ei.write(f"{E[0][i]:>12.6f}{E_2_x[0][i]:>12.6f}{E_2_y[0][i]:>12.6f}{E_2_z[0][i]:>12.6f}{E_2_xy[0][i]:>12.6f}{E_2_yz[0][i]:>12.6f}{E_2_zx[0][i]:>12.6f}")


        # # Now calcualting the absorption
        # abs_x = []
        # abs_y = []
        # abs_z = []
        # abs_E = []
        # with open(f"{OPTICS_DIR}/ABSORPTION.dat", "r+") as abs:
        #     for i,line in enumerate(abs):
        #         if i==0:continue
        #         split_line = line.split()
        #         if len(split_line) != 7: continue
        #         abs_E.append(float(split_line[0]))
        #         abs_x.append(float(split_line[1]))
        #         abs_y.append(float(split_line[2]))
        #         abs_z.append(float(split_line[3]))

        # print(abs_x[3])
        # plt.figure()
        # plt.plot(abs_E, abs_x, label=f"$XX$")
        # plt.plot(abs_E, abs_y, label=f"$YY$")
        # plt.plot(abs_E, abs_z, label=f"$ZZ$")
        # plt.legend()
        # # plt.title(f"Absorption vs Energy")
        # plt.xlabel("Energy (eV)")
        # plt.ylabel("Absorption ($cm^{-1})$")
        # plt.xlim([0, 3])
        # plt.ylim([0, 1000000])
        # plt.savefig(f"absvsE.png")
        # plt.show()
        

