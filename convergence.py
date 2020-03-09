import os
import datetime
import sys
import subprocess
import time
import matplotlib
# matplotlib.use('Agg')  # no UI backend required if working in the wsl without a UIs
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np
from ebk import eVA32GPa

class LatticeConstantOptimize():
    def __init__(self, a, E, name="a_optimization"):
        """
        Takes the lattice constants[list in Angstroms], energies[list in eV] and a name[for .out files] as inputs
        """
        self.a = a
        self.e = E  # everything is in eV now
        self.v = [0.25*n**3 for n in self.a]   # There are four primitive cells in a single conventional cell this is for a diamond/zinc blende structure
        self.name = name
        #Any other parameters you can set here:
        self.graph_title = "Optimization of lattice constant"

    def compute_fit(self):
        """
        This method calculates the fit parameters, optimized lattice costant/ volume / Bulk modulus.
        """
        self.z = np.polyfit(self.a, self.e, 2)  # Getting the fit parameters
        self.f = np.poly1d(self.z)  ## Getting the new function
        self.x_fit = np.linspace(self.a[0], self.a[-1], 100)
        self.y_fit = self.f(self.x_fit)

        # Similarly for the volume
        self.vz = np.polyfit(self.v, self.e, 2)  # Getting the fit parameters
        self.vf = np.poly1d(self.vz)  ## Getting the new function
        self.v_x_fit = np.linspace(self.v[0], self.v[-1], 100)
        self.v_y_fit = self.vf(self.v_x_fit)

        # Getting the minimum energy
        self.E_optimized = min(self.y_fit)
        self.E_optimized_printable = self.E_optimized.astype(np.float)

        # Getting the optimized lattice constant
        self.min_index = np.argmin(self.y_fit)
        self.a0_optimized = self.x_fit.flat[self.min_index]
        self.v0_optimized = self.v_x_fit.flat[self.min_index]      # There are four primitive cells in a single conventional cell

        # Calculations
        # Getting the double derivative using a 2nd degree polynomial
        self.dda0 = 2*self.z[0]#.flat[0]
        self.ddv0 = 2*self.vz[0]#.flat[0]
        self.B = eVA32GPa(self.v0_optimized*self.ddv0)  # 1 eV/Angstrom3 = 160.21766208 GPa

    def plot(self):
        plt.rcParams["figure.figsize"] = (14,9)
        fit_label = f"a$_0$: {round(self.a0_optimized,3)} $\\AA$, B$_0$: {round(self.B,3)} GPa, $\\Omega_0$:{round(self.v0_optimized,3)} $\\AA^3$"
        plt.ylabel('Total Energy (eV)')
        plt.xlabel('Lattice Constant ($\\AA$)')
        plt.plot(self.a, self.e, 'x', label="Data")
        plt.plot(self.x_fit, self.y_fit, '--', label=fit_label)
        plt.title(self.graph_title)
        plt.legend(loc='best')
        plt.savefig(f"{self.name}.pdf")
        plt.show()

class E_cut_Optimize():
    def __init__(self, E_cut, Energies, labels, name="", num_of_atoms = 2):
        """
        This initializes the E_cut optimization
        Inputs
        E_cut  # List of lists
        Energies  # List of lists
        label: List of strings
        """
        self.cut_off = E_cut  # This should be an array of arrays
        self.final_energies = []
        self.labels = labels
        for E in Energies:
            self.final_energies.append(np.array([x/num_of_atoms for x in E])) # Converting to per energies per atom
        self.name = name  # This part will be added to the file name
        self.final_DEs = []

        for E in self.final_energies:
            E0 = np.amin(E)
            self.final_DEs.append([(E1 - E0)*1000 for E1 in E])  # Converting in to meV
        self.graph_title = ""

    def plot(self, diff = True, MP = False, R = False):
        plt.rcParams["figure.figsize"] = (14,9)
        if diff == True:
            # self.final_DEs = [x*1000 for x in self.final_DEs]   # Converting in to meV
            E_to_plot = self.final_DEs
            # plt.plot(self.cut_off, self.final_DEs, 'x-')
            plt.ylabel("$\Delta$ E (meV/atom)")
        else:
            E_to_plot = self.final_energies
            # plt.plot(self.cut_off, self.final_energies, 'x-')
            plt.ylabel("Total Energy (eV/atom)")

        for x in range(0, len(self.cut_off)):
            plt.plot(self.cut_off[x], E_to_plot[x], 'x-', label = self.labels[x])

        plt.legend()
        # plt.xlim(28, 40)
        # plt.ylim(-3176.5, -3175.5)
        if MP == True:
            plt.title(f"Convergence SCF for wavefunction K-Grid cutoff {self.graph_title}")
            plt.xlabel("Monkhorst-Pack grid (3D)")
            plt.savefig(f"SCFconvergence_{self.name}_Diff{diff}_MPGrid.pdf")
        elif R ==  True:
            plt.title(f"Convergence SCF for $\rho$ cutoff {self.graph_title}")
            plt.xlabel("$\rho$ cutoff (Kinetic energy cutoff for charge density) (Ry)")
            plt.savefig(f"SCFconvergence_{self.name}_Diff{diff}_Rho.pdf")
        else:
            plt.title(f"Convergence SCF for wavefunction Kinetic Energy cutoff {self.graph_title}")
            plt.xlabel("Wave function cutoff (Ry)")
            plt.savefig(f"SCFconvergence_{self.name}_Diff{diff}_KE.pdf")
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

