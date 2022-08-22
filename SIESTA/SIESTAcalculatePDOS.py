"""
Reads both .DOS and .PDOS files. Can calculate band gaps. use the Read_PDOS class and load() functions to load() plots
Has the ability to pin the graphs at the homo/lumo or fermi level.
"""

import enum
import os
import matplotlib
import sys
import subprocess
import pathlib
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np
from matplotlib import gridspec
from ebk.SIESTA.SIESTAOutFileReader import SiestaReadOut

class Read_PDOS():
    def __init__(self, figure_name = "PDOS"):
        """
        self.pin = Can be homo/lumo/Ef or vacuum: homo and lumo work differently and pins the "lead" at zero with all others following. We can have multiple leads by setting in the load method and anything that follows will move respective to the last lead
        """

        self.system_labels = []
        self.file_names = []
        self.fermi_levels = []
        self.vacuum = []
        self.figure_name = figure_name
        self.Species_per_file = []
        self.show_band_gaps = True
        self.x = []
        self.y_up = []
        self.y_dn = []
        self.Egs = []
        self.set_y_range = False
        self.show_ef_as_line = False
        self.ylim_low = 0
        self.ylim_high = 50
        self.set_x_range = True
        self.xlim_low = -3
        self.xlim_high = 3
        self.orbital_labels = []
        self.plt_title = f""
        self.plt_ylabel = "PDoS"
        self.kwargs = []  # for all the matplotlib pltotting kwargs
        self.plt_show = True
        self.onsetup_xval = []
        self.onsetdown_xval = []
        self.onset_threshold = 0
        self.HOMOs_xval = []
        self.LUMOs_xval = []
        self.Lead_HOMOs_xval = []
        self.Lead_LUMOs_xval = []
        self.HOMOs_onset_xval = []
        self.LUMOs_onset_xval = []
        self.Lead_HOMOs_onset_xval = []
        self.Lead_LUMOs_onset_xval = []
        self.HOMO_Energies = []
        self.LUMO_Energies = []
        self.HOMO_from_vac = []
        self.pin = "vacuum"      # Can be homo/lumo/Ef or vacuum: homo and lumo work differently and pins the "lead" at zero with all others following. We can have multiple leads by setting in the load method and anything that follows will move respective to the last lead

    def process(self, system_label):
        """
        This functions automates loading the pdos file.
        system_label: String: System Label as string
        """
        self.system_labels.append(system_label)
        # Sn_out = SiestaReadOut("Sn")

        with open(f"{system_label}.PDOS", 'r') as file:
            for line in file:
                if "fermi_energy" in line:
                    self.fermi_levels.append(float(line.split("<")[-2].split(">")[-1]))

    def get_band_gaps(self):
        return self.Egs

    def load(self, file_name, Ef, label, lead = False, vacuum = 0, **kwargs):
        """
        This method will load the file and also calculate band gap when loaded. 
        Don't have to do the processing if using load (this method)
        lead = T/F Any energy shifts to move homo/lumo levels will follow the value at the lead.
                If no leads are set then first line will be the lead
                self.pin = "vacuum"      # Can be homo/lumo/Ef or vacuum: homo and lumo work differently and pins the "lead" 
                at zero with all others following. We can have multiple leads by setting in the load method and anything 
                that follows will move respective to the last lead

        **kwargs: ( matpotlib kwargs for lines and markers and the sort)
        """
        self.file_names.append(file_name)
        Ef = Ef
        self.orbital_labels.append(label)
        self.fermi_levels.append(float(Ef))
        self.kwargs.append(kwargs)
        tempx = []
        tempyup = []
        tempydn = []
        n = 0
        l = -1
        m = 9
        temp_n = 0
        temp_l = -1
        temp_m = 9
        self.vacuum.append(vacuum)

        with open(file_name, 'r') as file:
            for line in file:
                if "partial DOS for atom species" in line:
                    self.Species_per_file.append(line.strip().split(":")[-1].strip())

                if "#" not in line:
                    data = line.strip().split()
                    tempx.append(float(data[0]))
                    tempyup.append(float(data[1]))
                    if len(line) == 3:
                        tempydn.append(float(data[2]))
        self.x.append(tempx)
        self.y_up.append(tempyup)

        #finding the band gap
        diffs = [abs(x - Ef) for x in tempx]
        Ef_index = diffs.index(min(diffs))
        temp_negx = Ef_index
        temp_posx = Ef_index
        try:
            if tempyup[Ef_index] == 0:
                while tempyup[temp_negx] == 0:
                    temp_negx-=1
                while tempyup[temp_posx] == 0:
                    temp_posx+=1
                Eg = tempx[temp_posx-1] - tempx[temp_negx+1]
            self.HOMOs_xval.append(tempx[temp_negx+1])
            self.LUMOs_xval.append(tempx[temp_posx-1])
            if lead or self.Lead_HOMOs_xval == []: self.Lead_HOMOs_xval.append(tempx[temp_negx+1])
            else: self.Lead_HOMOs_xval.append(self.Lead_HOMOs_xval[-1])
            if lead or self.Lead_LUMOs_xval == []: self.Lead_LUMOs_xval.append(tempx[temp_posx-1])
            else: self.Lead_LUMOs_xval.append(self.Lead_LUMOs_xval[-1])
        except:
            #Most likely there was no band gap
            self.HOMOs_xval.append(0)
            self.LUMOs_xval.append(0)
            if lead or self.Lead_HOMOs_xval == []: self.Lead_HOMOs_xval.append(0)
            else: self.Lead_HOMOs_xval.append(self.Lead_HOMOs_xval[-1])
            if lead or self.Lead_LUMOs_xval == []: self.Lead_LUMOs_xval.append(0)
            else: self.Lead_LUMOs_xval.append(self.Lead_LUMOs_xval[-1])
            

        # This is for the calculations where you need to set a onset-threshold other than 0 to see where the DOS starts
        temp_negx = Ef_index
        temp_posx = Ef_index
        Eg = 0
        try:
            #We use a try catch since there might not be a band gap (conducting)
            if tempyup[Ef_index] == 0:
                while tempyup[temp_negx] <= self.onset_threshold:
                    temp_negx-=1
                while tempyup[temp_posx] <= self.onset_threshold :
                    temp_posx+=1
                Eg = tempx[temp_posx-1] - tempx[temp_negx+1]
            self.HOMOs_onset_xval.append(tempx[temp_negx+1])
            self.LUMOs_onset_xval.append(tempx[temp_posx-1])
            if lead or self.Lead_HOMOs_onset_xval == []: self.Lead_HOMOs_onset_xval.append(tempx[temp_negx+1])
            else: self.Lead_HOMOs_onset_xval.append(self.Lead_HOMOs_onset_xval[-1])
            if lead or self.Lead_LUMOs_onset_xval == []: self.Lead_LUMOs_onset_xval.append(tempx[temp_posx-1])
            else: self.Lead_LUMOs_onset_xval.append(self.Lead_LUMOs_onset_xval[-1])
            self.Egs.append(Eg)
        except:
            #Most likely there was no band gap
            self.HOMOs_onset_xval.append(0)
            self.LUMOs_onset_xval.append(0)
            if lead or self.Lead_HOMOs_onset_xval == []: self.Lead_HOMOs_onset_xval.append(0)
            else: self.Lead_HOMOs_onset_xval.append(self.Lead_HOMOs_onset_xval[-1])
            if lead or self.Lead_LUMOs_onset_xval == []: self.Lead_LUMOs_onset_xval.append(0)
            else: self.Lead_LUMOs_onset_xval.append(self.Lead_LUMOs_onset_xval[-1])
            self.Egs.append(Eg)

        self.HOMO_Energies.append(tempx[temp_negx+1])
        self.LUMO_Energies.append(tempx[temp_posx-1])

    def prepare(self):
        """Here we condition the data"""
        if self.pin.lower() == "homo":
            for i,line in enumerate(self.x):
                self.x[i] = [x - self.Lead_HOMOs_xval[i] for x in line]
        elif self.pin.lower() == "lumo":
            for i,line in enumerate(self.x):
                self.x[i] = [x - self.Lead_LUMOs_xval[i] for x in line]
        elif "vac" in self.pin.lower():
            for i,line in enumerate(self.x):
                self.x[i] = [x - self.vacuum[i] for x in line]
        elif "ef" in self.pin.lower():
            for i,line in enumerate(self.x):
                self.x[i] = [x - self.fermi_levels[i] for x in line]

    def plot(self):
        self.prepare()
        plt.figure()
        for i in range(0,len(self.x)):
            try:
                if self.show_band_gaps:
                    plt.plot(self.x[i], self.y_up[i], linewidth=1, label=f"{self.orbital_labels[i]:<6}\t$E_g$ = {self.Egs[i]:5> 2.3f} eV", **self.kwargs[i])
                else:
                    plt.plot(self.x[i], self.y_up[i], linewidth=1, label=f"{self.orbital_labels[i]:<6}", **self.kwargs[i])
            except:
                # If there are mutiple species or nothing picked up.
                if self.show_band_gaps:
                    plt.plot(self.x[i], self.y_up[i], linewidth=1, label=f"{self.orbital_labels[i]:<6}\t$E_g$ = {self.Egs[i]:> 2.3f} eV", **self.kwargs[i])
                else:
                    plt.plot(self.x[i], self.y_up[i], linewidth=1, label=f"{self.orbital_labels[i]:<6}", **self.kwargs[i])
            if self.show_ef_as_line:
                plt.axvline(self.fermi_levels[i])
        

        if self.pin.lower() == "ef": plt.xlabel(f'Energy in eV (E - E$_f$)')
        elif self.pin.lower() == "homo":  plt.xlabel(f'Energy in eV (E - E$_v$)')
        elif self.pin.lower() == "lumo": plt.xlabel(f'Energy in eV (E - E$_c$)')
        elif "vac" in self.pin.lower() : plt.xlabel(r'Energy in eV (E - E$_{vac}$)')

        plt.text(-3, 100, r'(b)', fontsize=12)
        plt.ylabel(f'{self.plt_ylabel}')
        plt.legend(loc='upper right')
        # plt.legend()
        # plt.title(f"{self.plt_title}")
        # plt.xticks(np.arange(-3,3, step=0.5))  # Set label locations.

        if self.set_x_range == True:
            plt.xlim([self.xlim_low,self.xlim_high])
        if self.set_y_range == True:
            plt.ylim([self.ylim_low,self.ylim_high])
        plt.savefig(f"{self.figure_name}.pdf")
        plt.tight_layout()
        if self.plt_show: plt.show()
        plt.close()

    def calculate_offsets(self):
        vals = []
        print("HOMO Offsets calculated from onset_threshold:",self.onset_threshold)
        for i in range(0,len(self.Lead_HOMOs_onset_xval)):
            vals.append(self.HOMOs_onset_xval[i] - self.Lead_HOMOs_onset_xval[i])
        zipped = zip(self.orbital_labels, vals)
        res = sorted(zipped, key = lambda x: x[1])

        for i,v in enumerate(vals):
            print(f"Shift: {self.orbital_labels[i]}: {vals[i]}")
        print(f"Ordered printing for ease of use")
        for i,v in enumerate(vals):
            print("Shift:", res[i])
        return vals

    def get_vacuum_subtracted_homo(self):
        """
        Returns the vacuum subtracted HOMO values and also prints them out
        """
        h_v = []
        for i,v in enumerate(self.HOMO_Energies):
            h_v.append(v - self.vacuum[i])
            print(f"HOMO-vac: {self.orbital_labels[i]}: {h_v[-1]}")
        return h_v