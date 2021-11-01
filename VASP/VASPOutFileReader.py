"""
This file reads the the file system_label.out and extracts/calculates data/values from it
"""
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
                    else: # This is because there might be partial occupancies
                        if E<= self.lowest_conduction[1]:
                            self.lowest_conduction[1] = E
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

