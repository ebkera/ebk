import os
import datetime
import sys
import subprocess
import time
import matplotlib
matplotlib.use('Agg')  # no UI backend required if working in the wsl without a UI
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np

class LatticeConstantOptimize():
    def __init__(self, a, E, name="a_optimization"):
        """
        Takes the lattice constants[list in Angstroms], energies[list in eV] and a name[for .out files] as inputs
        """
        self.a = a
        self.e = E
        self.name = name
                
        #Any other parameters you can set here
        self.graph_title = "Optimization of lattice constant"

    def compute_fit(self):
        self.z = np.polyfit(self.a, self.e, 3)  # Getting the fit parameters
        self.f = np.poly1d(self.z)  ## Getting the new function
        self.x_fit = np.linspace(self.a[0], self.a[-1], 100)
        self.y_fit = self.f(self.x_fit)

        # Getting the minimum energy
        self.E_optimized = min(self.y_fit)
        self.E_optimized_printable = self.E_optimized.astype(np.float)

        # Getting the optimized lattice constant
        self.min_index = np.argmin(self.y_fit)
        self.a0_optimized = self.x_fit.flat[self.min_index]

    def plot(self):
        fit_label = "Fit (" + str(round(self.a0_optimized,3)) + " $\\AA$)"
        plt.ylabel('Total Energy (eV)')
        plt.xlabel('Bond length ($\\AA$)')
        plt.plot(self.a, self.e, 'x', label="Data")
        plt.plot(self.x_fit, self.y_fit, '--', label=fit_label)
        plt.title(self.graph_title)
        plt.legend(loc='best')
        plt.savefig(f"{self.name}.pdf")

class E_cut_Optimize():
    def __init__(self, E_cut, E, name="", per_atoms = 2):
        self.cut_off = E_cut
        self.final_energies = [x/per_atoms for x in E]
        self.name = name
        self.final_energies = np.array(self.final_energies)*13.6056980659 # Converting to eV
        E0 = np.amin(self.final_energies)
        self.final_DEs = [E - E0 for E in self.final_energies]
        self.graph_title = ""

    def plot(self, diff = True, MP = False):
        plt.rcParams["figure.figsize"] = (14,9)
        if diff == True:
            plt.plot(self.cut_off, self.final_DEs, 'x-')
            plt.ylabel("$\Delta$ E (eV/atom)")
        else:
            plt.plot(self.cut_off, self.final_energies, 'x-')
            plt.ylabel("Total Energy (eV/atom)")
        # plt.xlim(28, 40)
        # plt.ylim(-3176.5, -3175.5)
        if MP == True:
            plt.title(f"Convergence SCF for wavefunction K-Grid cutoff {self.graph_title}")

        else:
            plt.title(f"Convergence SCF for wavefunction Kinetic Energy cutoff {self.graph_title}")
        plt.xlabel("Wave function cutoff (Ry)")
        plt.savefig(f"SCFconvergence_{self.name}_Diff{diff}_MP{MP}.pdf")
        plt.show()

# class K_cut_Optimize():
#     def __init__(self, k_cut, E, name="MP_optimization", per_atoms = 2):
#         self.MP_kpoints = k_cut
#         self.final_energies_kp = [x/per_atoms for x in E]
#         self.name = name
#         self.final_energies_kp = np.array(self.final_energies_kp)*13.6056980659 # Converting to eV
#         E0 = np.amin(self.final_energies_kp)
#         self.final_DEs = [E - E0 for E in self.final_energies_kp]

#     def plot(self, diff = True):
#         if diff == True:
#             plt.plot(self.MP_kpoints, self.final_DEs, 'x-')
#             plt.ylabel("$\Delta$ E (eV)")
#         else:
#             plt.plot(self.MP_kpoints, self.final_energies, 'x-')
#             plt.ylabel("Total Energy (eV)")
#         # plt.xlim(28, 40)
#         # plt.ylim(-3176.5, -3175.5)
#         plt.title("Convergence SCF for MP grid")
#         plt.xlabel("Wave function cutoff in Ry")
#         plt.savefig(f"Kconvergence_{self.name}.pdf")
#         plt.show()

if __name__ == "__main__":
    x = [55, 60, 65, 75, 80, 85]
    E = [-233.44606483, -233.44617876, -233.44627537, -233.44633003, -233.44633229, -233.44633816]
    O = E_cut_Optimize(x,E)
    O.graph_title = "For MP 15"
    # O.plot(False)
    O.plot(True)

    f55_10_300 = -233.44602856
    f60_10_300 = -233.44614256
    f65_10_300 = -233.44623909
    f55_15_300 = -233.44606483
    f60_15_300 = -233.44617876
    f65_15_300 = -233.44627537
    f75_15_300 = -233.44633003
    f80_15_300 = -233.44633229
    f85_15_300 = -233.44633816
    f65_20_300 = -233.44921409

    kcut = E_cut_Optimize([15, 20],[f65_15_300, f65_20_300])
    O.graph_title = "For KE cut 65"
    # # kcut.plot(True, True)
    # kcut.plot(False, True)
