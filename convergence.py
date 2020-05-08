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
    def __init__(self, E_cut, Energies, num_of_atoms, label, name=""):
        """
        This initializes the E_cut optimization
        Inputs
        E_cut (list): List of cutoff values can be K cut or KE_cut
        Energies (list): List of Energy values
        label: List of strings
        num_of_atoms: (integer) The number of atoms that the molecule/primitive unit cell has
        """
        self.cut_off = [E_cut]  # This should be an array of arrays
        self.Energies = [Energies]
        self.num_of_atoms = [num_of_atoms]
        self.labels = [label]
        self.per_atom = True
        self.name = name  # This part will be added to the file name
        self.graph_title = "SCF Convergence"
        self.showplot = True
        self.xlabel = ""
        self.ylabel = ""

    def add_to_plot(self, E_cut, Energies, num_of_atoms, label):
        """This function adds to extra plots to the same graph. See __init__ doc sting for varable references. They are same."""
        self.cut_off.append(E_cut)
        self.Energies.append(Energies)
        self.num_of_atoms.append(num_of_atoms)
        self.labels.append(label)

    def plot(self, diff = True, MP = False, R = False):
        """
        |This is the function that plots everything.
        |By default plots the Basis wavefucniton Kinetic energy convergence
        |diff: (bool) If true differences in energy is plotted
        |MP  : (bool) If True plots the convergence of Monkhorst-Pack grid
        |R   : (bool) If True plots the convergence of the charge density cutoff
        """
        self.final_energies = []
        self.final_DEs = []
        self.R = R
        for E in self.Energies:
            index = self.Energies.index(E)
            if self.per_atom:
                self.final_energies.append(np.array([x/self.num_of_atoms[index] for x in E])) # Converting to energies per atom
            else :
                self.final_energies.append(np.array([x for x in E])) # Not converting Converting to per energies per atom

        for E in self.final_energies:
            E0 = np.amin(E)
            self.final_DEs.append([(E1 - E0)*1000 for E1 in E])  # Converting in to meV

        plt.rcParams["figure.figsize"] = (14,9)
        plt.figure()

        # Setting the y label according to all the variables
        if diff == True:
            E_to_plot = self.final_DEs
            if self.per_atom:
                plt.ylabel(f"$\Delta$ E (meV/atom) {self.ylabel}")
            else:
                plt.ylabel(f"$\Delta$ E (meV) {self.ylabel}")
        else:
            E_to_plot = self.final_energies
            if self.per_atom:
                plt.ylabel(f"Total Energy (eV/atom) {self.ylabel}")
            else:
                plt.ylabel(f"Total Energy (eV) {self.ylabel}")

        for x in range(0, len(self.cut_off)):
            plt.plot(self.cut_off[x], E_to_plot[x], 'x-', label = self.labels[x])

        plt.legend()
        if MP == True:
            plt.title(f"{self.graph_title}")
            plt.xlabel(f"Monkhorst-Pack grid (3D) {self.xlabel}")
            plt.savefig(f"SCFconvergence_{self.name}_Diff{diff}_MPGrid_peratom{self.per_atom}.pdf")
        elif self.R == True:
            plt.title(f"{self.graph_title}")
            plt.xlabel(f"$\\rho$ cutoff (Kinetic energy cutoff for charge density) (Ry) {self.xlabel}")
            plt.savefig(f"SCFconvergence_{self.name}_Diff{diff}_Rho_peratom{self.per_atom}.pdf")
        else:
            plt.title(f"{self.graph_title}")
            plt.xlabel(f"Wave function cutoff (Ry) {self.xlabel}")
            plt.savefig(f"SCFconvergence_{self.name}_Diff{diff}_KE_peratom{self.per_atom}.pdf")
        if self.showplot: plt.show()

if __name__ == "__main__":
    pass

