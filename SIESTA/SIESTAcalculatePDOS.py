import os
import matplotlib
# matplotlib.use('Agg')  # no UI backend required if working in the wsl without a UI
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
        self.system_labels = []
        self.file_names = []
        self.fermi_levels = []
        self.figure_name = figure_name
        self.Species_per_file = []
        self.x = []
        self.y_up = []
        self.y_dn = []
        self.set_y_range = False
        self.ylim_low = 0
        self.ylim_high = 50
        self.set_x_range = True
        self.xlim_low = -3
        self.xlim_high = 3
        self.orbital_labels = []
        self.plt_title = f"{self.figure_name}"
        self.plt_ylabel = "PDOS"
        self.kwargs = []  # for all the matplotlib pltotting kwargs

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

    def load(self, file_name, Ef, label, **kwargs):
        """
        **kwargs: ( matpotlib kwargs for lines and markers and the sort)
        """
        self.file_names.append(file_name)
        Ef = Ef
        tempx = []
        tempyup = []
        tempydn = []
        n = 0
        l = -1
        m = 9
        temp_n = 0
        temp_l = -1
        temp_m = 9

        with open(file_name, 'r') as file:
            for line in file:
                if "partial DOS for atom species" in line:
                    self.Species_per_file.append(line.strip().split(":")[-1].strip())
                # if "Add data for atom_index" in line:
                #     new_m = line.strip().split()[-1].strip()
                #     new_l = line.strip().split()[-2].strip()
                #     new_n = line.strip().split()[-3].strip()

                #     if temp_n == 0: temp_n = new_n
                #     elif temp_n != int(new_n):

                if "#" not in line:
                    data = line.strip().split()
                    # print(line.strip().split()) X_val = 
                    # print(line)
                    tempx.append(float(data[0])-float(Ef))
                    tempyup.append(float(data[1]))
                    if len(line) == 3:
                        tempydn.append(float(data[2]))
        self.x.append(tempx)
        self.y_up.append(tempyup)
        self.orbital_labels.append(label)
        self.fermi_levels.append(float(Ef))
        self.kwargs.append(kwargs)

    def prepare(self):
        pass

    def plot(self):
        plt.figure()
        for i in range(0,len(self.x)):
            if "BDT" in self.orbital_labels[i]:
                try:
                    plt.plot(self.x[i], self.y_up[i], "--", linewidth=1, label=f"{self.orbital_labels[i]}", **self.kwargs[i])
                except:
                    plt.plot(self.x[i], self.y_up[i], "--", linewidth=1, label=f"{self.orbital_labels[i]}", **self.kwargs[i])
            else:
                try:
                    plt.plot(self.x[i], self.y_up[i], linewidth=1, label=f"{self.orbital_labels[i]}", **self.kwargs[i])
                except:
                    # If there are mutiple species or nothing picked up.
                    plt.plot(self.x[i], self.y_up[i], linewidth=1, label=f"{self.orbital_labels[i]}", **self.kwargs[i])
        plt.xlabel('Energy in eV (E - E$_f$)')
        # plt.text(-3, 100, r'(b)', fontsize=12)
        plt.ylabel(f'{self.plt_ylabel}')
        plt.legend(loc='upper right')
        # plt.legend()
        plt.title(f"{self.plt_title}")
        if self.set_x_range == True:
            plt.xlim([self.xlim_low,self.xlim_high])
        if self.set_y_range == True:
            plt.ylim([self.ylim_low,self.ylim_high])
        plt.savefig(f"{self.figure_name}.pdf")
        plt.show()
        plt.close()

        # plt.figure()
        # for Set in args:
        #     if Set.band_gap != None:
        #         plt.plot(Set.band_data_x, Set.band_data_y, linewidth=0.5, label=f"{Set.name}, $E_g$: {Set.band_gap:.3f}, $E_c-E_f$: {Set.gap_high:.3f}, $E_f - E_v$: {-Set.gap_low:.3f} eV")
        #     else:
        #         plt.plot(Set.band_data_x, Set.band_data_y, linewidth=0.5, label=f"{Set.name}, $E_g$: None")
        # plt.yscale('log')
        # plt.xlabel('Energy in eV (E - E$_f$)')
        # plt.ylabel('DOS')
        # # plt.legend(loc='upper left')
        # plt.legend()
        # plt.title("DOS")
        # plt.savefig(f"Comparison_log_{figure_name}.pdf")
        # plt.close()