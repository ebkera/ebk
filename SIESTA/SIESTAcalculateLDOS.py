"""
Plot the LDOS from files when a folder is given.

Notes:
Variable/attribute names may be non-pythonic since we are also making the script compatible with SIESTA. Camel case if prefered.

Works as follows:
    Sets the folder and file name and also the SystemLabel
    reads in the positions of the atoms and the species uses ase of this
    Reads in the PDOS files that are availabale and reads off of the labeling sceme SystemLabel.1.PDOS
    Reads in as user input the number of bins in the xyz directions that are required.

    Input paramers:
    folder_path
    SystemLabel
    PlotAlongAxis
    NumberOfBins
"""


from ebk.SIESTA.SIESTAOutFileReader import SiestaReadOut
from ebk.SIESTA.SIESTAcalculatePDOS import Read_PDOS
from ase.io import read, write
from matplotlib.pyplot import axes, get
from ebk.MatMan import get_extreme_coordinates
from numpy.core.numeric import zeros_like
from math import log10, log
import numpy as np
import matplotlib
from matplotlib.colors import LogNorm
from matplotlib import colors

class LDOS():
    def __init__(self, *args, **kwargs):
        # Files required and their paths are set here
        self.folder_path = kwargs.get("folder_path")
        self.SystemLabel = kwargs.get("SystemLabel")
        self.PlotAlongAxis = kwargs.get("PlotAlongAxis", "z")
        self.NumberOfBins = kwargs.get("NumberOfBins", 100)
        self.plt_title = kwargs.get("PltTitle", f"LDOS: {self.SystemLabel}")
        self.plt_set_ylim = kwargs.get("SetYlim", False)
        self.plt_set_ylim_max = kwargs.get("SetYlim_max", 3)
        self.plt_set_ylim_min = kwargs.get("SetYlim_min", -3)
        self.log = kwargs.get("log", False)
        self.cmap = kwargs.get("cmap","coolwarm") #Other Options:"viridis","plasma","Pastel1","tab20"

        # Getting other parameters automatically by reading in files
        self.figure_name = f"{self.SystemLabel}"
        self.outfile = SiestaReadOut(f"{self.folder_path}/{self.SystemLabel}")
        self.atoms = read(f"{self.folder_path}/{self.SystemLabel}.STRUCT_OUT")
        self.NumberOfAtoms = len(self.atoms)
        self.Ef = self.outfile.get_fermi()
        self.Energy = []  # An empty first element to make the index start from 1
        self.AtomicPDOS = {}  # An empty first element to make the index start from 1

        # Getting the dimensions of the cell and atoms
        self.extreme_coordinates = get_extreme_coordinates(self.atoms)
        # print(self.extreme_coordinates)
        self.len_of_axes = [self.extreme_coordinates[0][1]-self.extreme_coordinates[0][0],self.extreme_coordinates[1][1]-self.extreme_coordinates[1][0],self.extreme_coordinates[2][1]-self.extreme_coordinates[2][0]]

        # Read in all the files that are Required for atomic PDOS
        # First setting up PDOS 
        all_pdos = Read_PDOS("All_PDOS")
        for atom_n in range(0,self.NumberOfAtoms):
            try:
                # file_name = f"{self.folder_path}/{self.SystemLabel}.{atom_n}.PDOS"
                file_name = f"{self.folder_path}/{self.SystemLabel}.PDOS.{atom_n}"
                all_pdos.load(file_name, self.outfile.get_fermi(), f"{self.SystemLabel}.{atom_n}")
                self.AtomicPDOS.update({atom_n:all_pdos.y_up[-1]})
                self.Energy.append(all_pdos.x[-1])
            except:
                pass

        # Should be modified to reflect the porper axis right now harcoded for axis z
        self.axis_len = self.len_of_axes[2]
        self.axis_extremes = [self.extreme_coordinates[2][0], self.extreme_coordinates[2][1]]
        self.bin_size = self.axis_len/self.NumberOfBins

        # selecting atoms in a bin
        bin_min = self.axis_extremes[0]
        bin_max = bin_min + self.bin_size
        self.bins = []
        for bin in range(0,self.NumberOfBins):
            pdos_for_bin = []
            for i,v in enumerate(self.atoms):
                if v.position[2] >= bin_min and v.position[2] < bin_max:
                    if pdos_for_bin == []:
                        try: pdos_for_bin = self.AtomicPDOS[i].copy()
                        except: pass # this means that there was no atoms fo that atomic number with a PDOS file
                    else:
                        for x in range(0,len(pdos_for_bin)):
                            try: pdos_for_bin[x] += self.AtomicPDOS[i][x]
                            except: pass
            
            if len(pdos_for_bin) != 3000:
                for x in range(0,3000):
                    pdos_for_bin.append(0)
            self.bins.append(pdos_for_bin)
            bin_min = bin_max
            bin_max = bin_max + self.bin_size
        
        lowest = 1000
        if self.log == True:
            throwaway = self.bins.copy()
            self.bins = []
            for y in throwaway:
                t = [0.000001 if x == 0 else x for x in y]
                self.bins.append(t)
            throwaway = self.bins.copy()
            self.bins = []
            for y in throwaway:
                t = [np.log(x) for x in y]
                minimum = min(t) 
                if minimum < lowest:
                    lowest = minimum
                self.bins.append(t)
            throwaway = self.bins.copy()
            self.bins = []
            for y in throwaway:
                t = [x-lowest for x in y]
                self.bins.append(t)

    def plot(self):
        # divide the length of the required axis into bins
        # group atoms into the relevant bins according to their position
        # add up the PDOS for all the atoms in a single bin
        # plot the bins (z axis lenth) vs E to give image
        #   smooth over bins so that the image is continous

        import numpy as np
        import matplotlib.pyplot as plt
        xlist = np.linspace(self.axis_extremes[0],self.axis_extremes[1], self.NumberOfBins)
        ylist = self.Energy[-1]
        ylist = np.array(ylist)
        # ylist = ylist.transpose()
        print(xlist.shape,ylist.shape)
        X, Y = np.meshgrid(xlist,ylist)
        Z = np.array(self.bins).transpose()
        print(X.shape, Y.shape, Z.shape)
        fig,ax=plt.subplots(1,1)
        cp = ax.contourf(X, Y, Z, 20, cmap=self.cmap, nlevels=200)

        # This is the fix for the white lines between contour levels
        for c in cp.collections:
            c.set_edgecolor("face")

        norm= matplotlib.colors.Normalize(vmin=cp.cvalues.min(), vmax=cp.cvalues.max())
        sm = plt.cm.ScalarMappable(norm=norm, cmap = cp.cmap)
        sm.set_array([])
        if self.log:fig.colorbar(sm, ticks=cp.levels, label="log(DoS)")
        else:fig.colorbar(sm, ticks=cp.levels, label="DoS")
        ax.set_title(self.plt_title)
        ax.set_xlabel('Position ($\AA$)')
        if self.plt_set_ylim:
            ax.set_ylim((self.plt_set_ylim_min, self.plt_set_ylim_max))
        ax.set_ylabel('E (eV)')
        plt.tight_layout()
        if self.log: log_text="_log"
        else: log_text = ""
        plt.savefig(f"{self.figure_name}{log_text}_LDOS.pdf")
        plt.savefig(f"{self.figure_name}{log_text}_LDOS.png")
        # plt.show() 
