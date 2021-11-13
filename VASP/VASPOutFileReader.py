"""
This file reads the the file system_label.out and extracts/calculates data/values from it
"""
from os import P_WAIT, kill
from ebk.BandPlotter import BandPlotter 

class VASPReadOut():
    def __init__(self, out_folder, *args, **kwargs):
        self.out_folder = out_folder
        self.bands = []
        self.kpoints = []
        trigger = False
        self.k_dist=[0]
        self.highest_valance = [[0,0,0], -500]    # Will contain [[kpathpoint],E]
        self.lowest_conduction = [[0,0,0], 500]  # Will contain [[kpathpoint],E]
        self.valance_band_label_index = 0
        self.conduction_band_label_index = 0
        self.hsp = []
        self.hss = []
        self.Eg = "notset"
        self.draw_band_edge_lines = False
        self.spin_orbit = None
        self.Ef = 0
        
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

                # print(line)  # Left here for debugging 
                if 'k-point' in line and 'plane waves' in line:
                    k = line.split()
                    k2 = [float(k[3]),float(k[4]),float(k[5])]
                    self.kpoints.append(k2)
                    # print(k2)

                if 'E-fermi' in line:
                    k = line.split()
                    E_fermi = float(k[2])
                    self.Ef = E_fermi
                    # E_fermi = 0
                    # print(E_fermi)

                if " k-point     1 " in line:
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
                    band_label = int(f[0])-1
                    E = float(f[1])-E_fermi
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
                            self.valance_band_label_index = band_label
                    else: # This is because there might be partial occupancies
                        if E<= self.lowest_conduction[1]:
                            self.lowest_conduction[1] = E
                            self.conduction_band_label_index = band_label
                            # print(E)

        for x in range(1,len(self.kpoints)):
            dist_v=[0, 0, 0]
            for i in range(3):
                # print (i, x)
                dist_v[i] = self.kpoints[x][i] - self.kpoints[x-1][i]
            self.k_dist.append((dist_v[0]**2+dist_v[1]**2+dist_v[2]**2) + self.k_dist[x-1])

        with open(f"{self.out_folder}/KPOINTS", "r+") as KPOINTS:
            hss_old = "empty"
            for line in KPOINTS:
                line_s = line.split()

                if len(line_s) == 1:
                    try:
                        self.k_point_density = int(line_s[0])
                    except: pass

                # This will check for actual K points
                try: 
                    float(line_s[0])
                except: 
                    continue

                if len(line_s) == 4:
                    hss = line_s[3]
                    if hss == hss_old: continue
                    hss_old = hss
                    self.hss.append(hss)

                    for i,v in enumerate(self.hss):
                        if v.lower() == "gamma":
                            self.hss[i] = "$\Gamma$"
            self.hsp.append(self.k_dist[0])

            for x in range(1,len(self.hss)):
                self.hsp.append(self.k_dist[self.k_point_density*x-1])

        self.Eg = self.lowest_conduction[1] - self.highest_valance[1]

    def get_band_gap(self):
        return self.Eg

    def get_fermi_energy(self):
        return self.Ef

    def get_band_structure(self, title = "Band Diagram", file_name = None):
        """
        Sets initial settings and returns a BandPlotter type object for further manipulation but ready to plot
        """
        if file_name == None:
            file_name = self.out_folder
        
        x = BandPlotter(self.k_dist, self.bands, self.hsp, self.hss)
        x.file_name = f"{file_name}"
        x.title = title
        x.set_y_range=True
        x.same_band_colour = True
        x.hlines = True
        x.vlines = True
        # x.ylim_high = 2.5
        # x.ylim_low = -2.5
        # x.ylim_high = 3
        # x.ylim_low = -15
        x.saveas_extension = "png"

        if self.draw_band_edge_lines:
            x.extra_horizontal_lines = [[self.lowest_conduction[1], f"Conduction band edge"]]
            x.extra_horizontal_lines.append([self.highest_valance[1], f"Valance band edge"])
        return x


    def get_curvature_at_k_point(self, k_point: list[float,float,float]=False, k_index:int=False, direction:str=False, sigma:int=20, high_symmetry_point:bool=True, show_plot:bool=False)-> float:
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

            electron_effective_mass_left = scipy.constants.hbar**2/(2*fit_para_left_conduction[0]*scipy.constants.m_e)  # The two on the denominator is from the 2 that comes from differentiating a 2 order polynomial
            electron_effective_mass_right = scipy.constants.hbar**2/(2*fit_para_right_conduction[0]*scipy.constants.m_e)  # The two on the denominator is from the 2 that comes from differentiating a 2 order polynomial
            hole_effective_mass_left = scipy.constants.hbar**2/(2*fit_para_left_valance[0]*scipy.constants.m_e)  # The two on the denominator is from the 2 that comes from differentiating a 2 order polynomial
            hole_effective_mass_right = scipy.constants.hbar**2/(2*fit_para_right_valance[0]*scipy.constants.m_e)  # The two on the denominator is from the 2 that comes from differentiating a 2 order polynomial

            return_dict = {}
            return_dict.update({"fit_para_conduction_l":fit_para_left_conduction})
            return_dict.update({"fit_para_conduction_r":fit_para_right_conduction})
            return_dict.update({"fit_para_valance_l":fit_para_left_valance})
            return_dict.update({"fit_para_valance_r":fit_para_left_valance})
            return_dict.update({"me_l":electron_effective_mass_left})
            return_dict.update({"me_r":electron_effective_mass_right})
            return_dict.update({"mh_l":hole_effective_mass_left})
            return_dict.update({"mh_r":hole_effective_mass_right})
            return_dict.update({"bottom_of_conduction_l": fit_para_left_conduction[2]})
            return_dict.update({"bottom_of_conduction_r":fit_para_right_conduction[2]})
            return_dict.update({"top_of_valance_l":fit_para_left_valance[2]})
            return_dict.update({"top_of_valance_r":fit_para_right_valance[2]})

            p_era_LDA_ES = np.poly1d(fit_para_left_conduction)
            fit_curve_left_conduction = p_era_LDA_ES(left_side_kpoint_dist)
            p_right_conduction = np.poly1d(fit_para_right_conduction)
            fit_curve_right_conduction = p_right_conduction(right_side_kpoint_dist)
            
            p_left_valance = np.poly1d(fit_para_left_valance)
            fit_curve_left_valance = p_left_valance(left_side_kpoint_dist)
            p_right_valance = np.poly1d(fit_para_right_valance)
            fit_curve_right_valance = p_right_valance(right_side_kpoint_dist)
            import matplotlib.pyplot as plt 
            plt.plot(left_side_kpoint_dist, fit_curve_left_conduction, label=f"Left Conduction fit")
            plt.plot(left_side_kpoint_dist, left_side_conduction_E, label=f"Left Conduction")
            plt.plot(right_side_kpoint_dist, fit_curve_right_conduction, label=f"Right Conduction fit")
            plt.plot(right_side_kpoint_dist, right_side_conduction_E, label=f"Right Conduction")
            plt.plot(left_side_kpoint_dist, fit_curve_left_valance, label=f"Left Valance fit")
            plt.plot(left_side_kpoint_dist, left_side_valance_E, label=f"Left Valance")
            plt.plot(right_side_kpoint_dist, fit_curve_right_valance, label=f"Right Valance fit")
            plt.plot(right_side_kpoint_dist, right_side_valance_E, label=f"Right Valance")
            plt.legend()
            if show_plot: plt.show()

            return return_dict

        if direction:
            pass
