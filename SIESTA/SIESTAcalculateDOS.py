"""
This Files needs the eig2dos file to work
"""
import os
import matplotlib
matplotlib.use('Agg')  # no UI backend required if working in the wsl without a UI
import sys
import subprocess
import pathlib
from pathlib import Path


import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np
from matplotlib import gridspec
# import ebk.SIESTA.eig2dos as eig2dos

def latexify(inputlist):
    """This method will convert a list of strings into it's latex form"""
    for index, item in enumerate(inputlist):
        if item == "Gamma":
            inputlist[index] = '$\Gamma$'
    return inputlist

def plotter(figure_name = "DOS", *args, **kwargs):
    plt.figure()
    for Set in args:
        plt.plot( Set.band_data_x, Set.band_data_y, linewidth=0.5, label=f"{Set.name}, $E_g$: {Set.band_gap:.3f}, $E_c-E_f$: {Set.gap_high:.3f}, $E_f - E_v$: {-Set.gap_low:.3f} eV")
    plt.xlabel('Energy in eV (E - E$_f$)')
    plt.ylabel('DOS')
    # plt.legend(loc='upper left')
    plt.legend()
    plt.title(f"DOS")
    plt.savefig(f"Comparison_{figure_name}.pdf")
    plt.close()

    plt.figure()
    for Set in args:
        plt.plot( Set.band_data_x, Set.band_data_y, linewidth=0.5, label=f"{Set.name}, $E_g$: {Set.band_gap:.3f}, $E_c-E_f$: {Set.gap_high:.3f}, $E_f - E_v$: {-Set.gap_low:.3f} eV")
    plt.yscale('log')
    plt.xlabel('Energy in eV (E - E$_f$)')
    plt.ylabel('DOS')
    # plt.legend(loc='upper left')
    plt.legend()
    plt.title("DOS")
    plt.savefig(f"Comparison_log_{figure_name}.pdf")
    plt.close()

def plotter_vertical(*args, **kwargs):
    plt.figure()
    plt.subplot(121)   
    for Set in args:
        plt.plot( Set.band_data_y, Set.band_data_x, linewidth=0.4, label=Set.name)
    plt.ylabel('Energy in eV (E - E$_f$)')
    plt.xlabel('DOS')
    plt.legend(loc='upper right')
    plt.title("DOS comparison")
    # plt.savefig("Comparison.pdf")

    plt.subplot(122)
    for Set in args:
        plt.plot( Set.band_data_y, Set.band_data_x, linewidth=0.4, label=Set.name)
    plt.xscale('log')
    # plt.ylabel('Energy in eV (E - E$_f$)')
    plt.xlabel('DOS (log)')
    plt.legend(loc='upper right')
    plt.title("DOS (log) comparison")
    plt.savefig("Comparison_vertical.pdf")

def plotter_horizontal_vertical(*args, **kwargs):
    plt.figure()
    plt.subplot(121)   
    for Set in args:
        plt.plot( Set.band_data_x, Set.band_data_y, linewidth=0.4, label=Set.name)
    plt.xlabel('Energy in eV (E - E$_f$)')
    plt.ylabel('DOS')
    plt.legend(loc='upper left')
    plt.title("DOS comparison")

    plt.subplot(122)
    for Set in args:
        plt.plot( Set.band_data_y, Set.band_data_x, linewidth=0.4, label=Set.name)
    plt.ylabel('Energy in eV (E - E$_f$)')
    plt.xlabel('DOS')
    plt.legend(loc = 'upper right')
    plt.title("DOS comparison")
    plt.savefig("Comparison_horizontal_vertical.pdf")

def plotter_grid(*args, **kwargs):
    fig = plt.figure(figsize=(8, 6)) 
    gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 
    ax0 = plt.subplot(gs[0])
    for Set in args:
        ax0.plot( Set.band_data_x, Set.band_data_y, linewidth=0.4, label=Set.name) 
    # ax0.xlabel('Energy in eV (E - E$_f$)')
    # ax0.ylabel('DOS')
    ax0.legend(loc = 'upper left')
    ax1 = plt.subplot(gs[1])
    for Set in args:
        ax1.plot( Set.band_data_y, Set.band_data_x, linewidth=0.4, label=Set.name)
    # ax1.ylabel('Energy in eV (E - E$_f$)')
    # ax1.xlabel('DOS')
    ax1.legend(loc = 'upper right')
    plt.tight_layout()
    plt.title("DOS comparison")
    plt.savefig('grid_figure.pdf')

class Band():
    """
    Main class that deals with DOS. Dont know why Ii named it bands but retaining it for legacy reasons.
    """
    def __init__(self, name):
        self.name = name
        self.file_name_eig = name + ".EIG"
        self.file_name_dos = name + ".DOS"
        self.band_data_x = []
        self.band_data_y = []
        self.Efshift = False
        self.band_number_min = 0
        self.band_number_max = 0
        self.E_min = 0
        self.E_max = 0
        self.E_f = 0
        self.kpath_k = []
        self.kpath_symbol = []

    def automate(self):
        self.modify_eig_file(0.05, 1000, -3, 3)
        self.make_DOS_file()
        self.load_from_dos()

    def modify_eig_file(self, eta, ne, emin, emax):
        """Opens the EIG file and then saves the modified version of it
        eta  - I'd call it precision, put it less than 1. 
        Ne - number of eigenvalues, I advise to put it large. 
        Emin and Emax - interval of energies for which DOS will be calculated (eV). """
        file = open(self.file_name_eig,"r")
        rows = []
        for line in file:
            rows.append(line)
        file.close()

        # This is to catch if the first row is not a float. If it is will proceed. If had been edited we get a warning and will still replace the old parameters with new ones
        try:
            self.E_f = float(rows[0].strip())
        except:
            print("Warning: The " + self.file_name_eig + " file seems to have been processed previously. Setting parameters again...")
            self.E_f = float(rows[0].split()[0].strip())

        rows[0] = str(self.E_f) + " " + str(eta) + " " + str(ne) + " " + str(emin) + " " + str(emax) +"\n"

        file2 = open(self.file_name_eig,"w+")
        for row in rows:
            file2.write(row)
        file2.close()
    
    def make_DOS_file(self):
        self.system_version = sys.platform
        # self.eig2dos_path = os.path.abspath(__file__)
        self.eig2dos_path = pathlib.Path(__file__).parent.absolute()
        # print(f"System detected: {self.system_version}")
        # print(f"The required eig2dos file is expected to be in: {self.eig2dos_path}")
        # filename = Path("source_data/text_files/raw_data.txt")
        # self.eig2dos_path = Path(self.eig2dos_path) / "eig2dos"
        # print(f"Os.path {Path(self.eig2dos_path)}")
        if self.system_version == "win32":
            # print(f"wsl {self.eig2dos_path}")
            # os.system(f"wsl {self.eig2dos_path} < {self.file_name_eig} > {self.file_name_dos}")
            os.system("wsl ./eig2dos  < " + self.file_name_eig + " > " + self.file_name_dos)
        else:
            os.system("eig2dos  < " + self.file_name_eig + " > " + self.file_name_dos)


    def calculate_band_gap(self, calculate_using_peaks = False):
        """calculates teh band gaps"""
        if calculate_using_peaks:
            pass
        else:
            E_of_non_zeros = [i for i,v in enumerate(self.band_data_x) if self.band_data_y[i] == 0]
            E_of_non_zeros.append(E_of_non_zeros[-1]+2)  # just in case there is only one badn gap and therefpre there are no extra dataponts where it all changes so that the change can be detected.
            #calculating the band gap
            gap_low_index = E_of_non_zeros[0]
            gap_low = self.band_data_x[gap_low_index]
            for i in range(1,len(E_of_non_zeros)):
                if E_of_non_zeros[i]-E_of_non_zeros[i-1] != 1:
                    gap_high_index = E_of_non_zeros[i-1]
                    gap_high = self.band_data_x[gap_high_index]
                    if gap_low <= 0 and gap_high >= 0:
                        self.band_gap = gap_high - gap_low
                        self.gap_high = gap_high
                        self.gap_low = gap_low
                        # print(f"Band Gap for {self.name}: {self.band_gap} eV")
                        # print(f"Lower  edge: {gap_low} eV")
                        # print(f"Higher edge: {gap_high} eV")
                        break
                    gap_low_index = E_of_non_zeros[i]
                    gap_low = self.band_data_x[gap_low_index]

    def load_from_dos(self):
        """Opens the files to read"""
        file = open(self.file_name_dos,"r")
        rows = []
        for line in file:
            rows.append(line)
        file.close()

        # Delete the first 11 rows Since it just contains meta data
        for x in range(0,11):
            del rows[0]
        
        #saving the data to x and y band lists
        for row in rows:
            data_point = row.split()
            self.band_data_x.append(float(data_point[0]))
            self.band_data_y.append(float(data_point[1]))

if __name__ == "__main__": 
    Sn29 = Band("Sn29")
    Sn29.automate()
    Sn47 = Band("Sn47")
    Sn47.automate()
    plotter(Sn29,Sn47)
    plotter_vertical(Sn29,Sn47)
    plotter_horizontal_vertical(Sn29,Sn47)
    plotter_grid(Sn29,Sn47)