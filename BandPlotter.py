import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from ase.dft.kpoints import resolve_custom_points, find_bandpath_kinks
from ase.dft.dos import DOS

# matplotlib.use('Agg')  # no UI backend required if working in the wsl without a UI


class BandPlotter():
    def __init__(self, x, y, k_locations, k_symbols):
        """
        |All the inputs needed for a band plot are set here
        |Inputs:
        |   x          : The k path distance as floats
        |   y          : This should be a list of lists with the lists being the individual bands
        |   k_locations: The locations of the High symmetry points as a list
        |   k_symbols  : The Latexified version of the K path High symmetry point symbols
        """
        self.file_name = "band_diagram"
        self.x = x
        self.y = y
        self.k_symbols = k_symbols
        self.k_locations = k_locations
        self.number_of_bands_to_plot = 30
        self.fermi_level = 0
        self.hlines = False
        self.vlines = False
        self.title = "Band diagram"
        self.set_y_range = False
        self.ylim_low = -5
        self.ylim_high = 5
        self.xlim_low = 0
        self.xlim_high = 0
        self.Ef_shift = 0
        self.y_to_plot = []
        self.same_band_colour = False
        self.band_colour = "b"
        self.new_fig = False

    def plot(self):
        """
        |All the features of the band plot are set here
        |Inputs: None
        """
        for i in range(0,len(self.y)):
            self.y_to_plot.append([x - self.Ef_shift for x in self.y[i]])
        
        # To get multiple band diagrams together just plot them on the same figure (with the same file name) and when ever you want a new figure with just the new plots
        #  just do new_fig = True. If you want to save individual figures just give a new file name (with new_fig = True).
        if self.new_fig == True:
            plt.figure()
        
        # Setting vertical lines
        if self.vlines == True:
            for xc in self.k_locations:  # Plotting a dotted line that will run vertical at high symmetry points on the Kpath
                plt.axvline(x=xc, linestyle='-', color='k', linewidth=0.1)

        # Setting y ranges
        if self.set_y_range == True:
            plt.ylim([self.ylim_low,self.ylim_high])

        # We plot the figure here
        for band in range(0,len(self.y)):
            if self.same_band_colour == True:
                plt.plot(self.x, self.y_to_plot[band], self.band_colour)
            else:
                plt.plot(self.x, self.y_to_plot[band])

        plt.xticks(self.k_locations, self.k_symbols)
        plt.xlabel("K path")
        plt.ylabel("Energy (eV)")
        plt.title(f"{self.title}")
        plt.savefig(f"Bands_{self.file_name}.pdf")
        plt.show()


class BandPlotterASE():
    # Under costruction!!
    def __init__(self, **kwargs):
        """
        If plotting multiple band plots here can do only same path.
        The idea of how to use this class is as below:
            1) you set the required parameters that you need to plot the graph
            2) You add all the bands you want to add by using the add_to_plot() method
            3) Plot all bands that you have added using above method by calling the plot() method
        """
        self.file_name = "band_diagram"
        self.number_of_bands_to_plot = 30
        self.fermi_level = 0
        self.hlines = False
        self.vlines = False
        self.title = "Band diagram"
        self.set_y_range = False
        self.ylim_low = -5
        self.ylim_high = 5
        self.xlim_low = 0
        self.xlim_high = 0
        self.Ef_shift = True
        self.y_to_plot = []
        self.x_to_plot = []
        self.labels = []
        self.same_band_colour = False
        self.band_colour = ["b", "g", "r", "c", "m", "y", "k"]
        self.new_fig = False
        self.k_locations = None
        self.k_symbols = None
        self.dots = kwargs.get("dots", False)

    def plot(self):
        """
        |All the features of the band plot are set here
        |Inputs: None
        """

        # Setting the dimensions of the saved image
        plt.rcParams["figure.figsize"] = (14,9)

        # Setting vertical lines
        if self.vlines == True:
            for xc in self.k_locations:  # Plotting a dotted line that will run vertical at high symmetry points on the Kpath
                plt.axvline(x=xc, linestyle='-', color='k', linewidth=0.1)

        # Setting horizontal line on the fermi energy
        if self.hlines == True:
            plt.axhline(linewidth=0.1, color='k')

        # Setting y ranges
        if self.set_y_range == True:
            plt.ylim([self.ylim_low, self.ylim_high])

        # We plot the figure here
        for structure in range(0,len(self.y_to_plot)):
            # print(f"structure: {}")
            for band in range(0,len(self.y_to_plot[structure])):
                if self.same_band_colour == True:
                    if self.labels[structure] != None:
                        plt.plot(self.x_to_plot[structure], self.y_to_plot[structure][band], self.band_colour[structure], label = self.labels[structure])
                        self.labels[structure] = None
                    else:
                        plt.plot(self.x_to_plot[structure], self.y_to_plot[structure][band], self.band_colour[structure], label = self.labels[structure])
                else:
                    plt.plot(self.x_to_plot[structure], self.y_to_plot[structure][band])


        # # This is only for testing with Pierre and should be removed once done
        # # We plot the figure here
        # for structure in range(0,len(self.y_to_plot)):
        #     # print(f"structure: {}")
        #     for band in range(0,len(self.y_to_plot[structure])):
        #         if self.same_band_colour == True:
        #             if self.labels[structure] != None:
        #                 plt.plot(self.x_to_plot[structure], self.y_to_plot[structure][band], ".", label = self.labels[structure])
        #                 self.labels[structure] = None
        #             else:
        #                 plt.plot(self.x_to_plot[structure], self.y_to_plot[structure][band], ".", label = self.labels[structure])
        #         else:
        #             plt.plot(self.x_to_plot[structure], self.y_to_plot[structure][band], ".")


        plt.xticks(self.k_locations, self.k_symbols)
        plt.xlabel("K path")
        plt.ylabel("Energy (eV)")
        plt.title(f"{self.title}")
        plt.legend()
        plt.savefig(f"Bands_{self.file_name}.pdf")
        plt.show()

    def add_to_plot(self, readoutfilesobj, label = None):
        """
        |Here you add individual plots that need to be plot and then just plot them with the plot() method
        |Use this method which is a part of the BandPlotterASE class you will have to give the bands to plot using a readoutfilesobj type of object
        """
        try:
            Ef = readoutfilesobj.atoms_objects[0].calc.get_fermi_level()
        except:
            print(f"add_to_plot: Could not read fermi level")
        try:
            kpts = readoutfilesobj.atoms_bands_objects[0].calc.get_ibz_k_points()
        except:
            print(f"add_to_plot: Could not read k points")

        # Test space for k path and k high symmetry points
        # print(kpts)
        path = readoutfilesobj.atoms_bands_objects[0].cell.bandpath(npoints=0)
        kinks = find_bandpath_kinks(readoutfilesobj.atoms_bands_objects[0].cell, kpts, eps=1e-5)  # These are the high symmetry points in use 
        pathspec = resolve_custom_points(kpts[kinks], path.special_points, eps=1e-5) # This gives the postions for the relevant high symmetry points
        path.kpts = kpts
        path.path = pathspec

        klengths = []
        for x in range(0, len(kpts)):
            if x == 0:
                kdist = np.sqrt((0.0 - kpts[x][0])**2 + (0.0 - kpts[x][1])**2 +(0.0 - kpts[x][2])**2)
                klengths.append(kdist)
            else:
                kdist = np.sqrt((kpts[x-1][0] - kpts[x][0])**2 + (kpts[x-1][1] - kpts[x][1])**2 + (kpts[x-1][2]- kpts[x][2])**2)
                klengths.append(kdist+klengths[x-1])
        self.k_locations = []

        for x in range(len(kinks)):
            self.k_locations.append(klengths[kinks[x]])

        self.k_symbols = []
        for x in pathspec:
            if x == "G":
                self.k_symbols.append("$\Gamma$")
            else:
                self.k_symbols.append(x)

        energies = []
        for s in range(readoutfilesobj.atoms_bands_objects[0].calc.get_number_of_spins()):
            energies.append([readoutfilesobj.atoms_bands_objects[0].calc.get_eigenvalues(kpt=k, spin=s) for k in range(len(kpts))])
        # print(f"lenght of kpoints: {range(len(kpts))}")
        Energy_to_plot = []
        if self.Ef_shift == True:
            for band in energies[0]:
                Energy_to_plot.append([E - Ef for E in band])
        tempMain = []
        temp = []
        for x in range(len(Energy_to_plot[0])):
            for kpoint in Energy_to_plot:
                temp.append(kpoint[x])
            tempMain.append(temp)
            temp = []
        self.y_to_plot.append(tempMain)
        self.x_to_plot.append(klengths)
        self.labels.append(label)

class DOSPlotterASE():
    # Under costruction!!
    def __init__(self, **kwargs):
        """
        If plotting multiple band plots here can do only same path.
        The idea of how to use this class is as below:
            1) you set the required parameters that you need to plot the graph
            2) You add all the bands you want to add by using the add_to_plot() method
            3) Plot all bands that you have added using above method by calling the plot() method
        """
        self.file_name = "DOS"
        self.fermi_level = 0
        self.hlines = False
        self.vlines = False
        self.title = "DOS"
        self.set_y_range = False
        self.ylim_low = -5
        self.ylim_high = 5
        self.xlim_low = 0
        self.xlim_high = 0
        self.Ef_shift = True
        self.y_to_plot = []
        self.x_to_plot = []
        self.labels = []
        self.same_band_colour = False
        self.band_colour = ["b", "g", "r", "c", "m", "y", "k"]
        self.new_fig = False
        self.k_locations = None
        self.k_symbols = None
        self.dots = kwargs.get("dots", False)

    def plot(self):
        """
        |All the features of the band plot are set here
        |Inputs: None
        """

        # Setting the dimensions of the saved image
        plt.rcParams["figure.figsize"] = (14,9)

        # Setting vertical lines
        if self.vlines == True:
            for xc in self.k_locations:  # Plotting a dotted line that will run vertical at high symmetry points on the Kpath
                plt.axvline(x=xc, linestyle='-', color='k', linewidth=0.1)

        # Setting horizontal line on the fermi energy
        if self.hlines == True:
            plt.axhline(linewidth=0.1, color='k')

        # Setting y ranges
        if self.set_y_range == True:
            plt.ylim([self.ylim_low, self.ylim_high])

        # We plot the figure here
        for i,v in enumerate(self.y_to_plot):
            plt.plot(self.x_to_plot[i], self.y_to_plot[i], self.band_colour[i], label = self.labels[i])

        # plt.xticks(self.k_locations, self.k_symbols)
        plt.ylabel("DOS")
        plt.xlabel("Energy (eV)")
        plt.title(f"{self.title}")
        plt.legend()
        plt.savefig(f"DOS_{self.file_name}.pdf")
        plt.show()

    def add_to_plot(self, readoutfilesobj, label = None):
        """
        |Here you add individual plots that need to be plot and then just plot them with the plot() method
        |Use this method which is a part of the BandPlotterASE class you will have to give the bands to plot using a readoutfilesobj type of object
        """
        # try:
        #     Ef = readoutfilesobj.atoms_objects[0].calc.get_fermi_level()
        # except:
        #     print(f"add_to_plot: Could not read fermi level")
        # try:
        #     kpts = readoutfilesobj.atoms_bands_objects[0].calc.get_ibz_k_points()
        # except:
        #     print(f"add_to_plot: Could not read k points")

        # # Test space for k path and k high symmetry points
        # # print(kpts)
        # path = readoutfilesobj.atoms_bands_objects[0].cell.bandpath(npoints=0)
        # kinks = find_bandpath_kinks(readoutfilesobj.atoms_bands_objects[0].cell, kpts, eps=1e-5)  # These are the high symmetry points in use 
        # pathspec = resolve_custom_points(kpts[kinks], path.special_points, eps=1e-5) # This gives the postions for the relevant high symmetry points
        # path.kpts = kpts
        # path.path = pathspec

        # klengths = []
        # for x in range(0, len(kpts)):
        #     if x == 0:
        #         kdist = np.sqrt((0.0 - kpts[x][0])**2 + (0.0 - kpts[x][1])**2 +(0.0 - kpts[x][2])**2)
        #         klengths.append(kdist)
        #     else:
        #         kdist = np.sqrt((kpts[x-1][0] - kpts[x][0])**2 + (kpts[x-1][1] - kpts[x][1])**2 + (kpts[x-1][2]- kpts[x][2])**2)
        #         klengths.append(kdist+klengths[x-1])
        # self.k_locations = []

        # for x in range(len(kinks)):
        #     self.k_locations.append(klengths[kinks[x]])

        # self.k_symbols = []
        # for x in pathspec:
        #     if x == "G":
        #         self.k_symbols.append("$\Gamma$")
        #     else:
        #         self.k_symbols.append(x)

        # energies = []
        # for s in range(readoutfilesobj.atoms_bands_objects[0].calc.get_number_of_spins()):
        #     energies.append([readoutfilesobj.atoms_bands_objects[0].calc.get_eigenvalues(kpt=k, spin=s) for k in range(len(kpts))])
        # # print(f"lenght of kpoints: {range(len(kpts))}")
        # Energy_to_plot = []
        # if self.Ef_shift == True:
        #     for band in energies[0]:
        #         Energy_to_plot.append([E - Ef for E in band])
        # tempMain = []
        # temp = []
        # for x in range(len(Energy_to_plot[0])):
        #     for kpoint in Energy_to_plot:
        #         temp.append(kpoint[x])
        #     tempMain.append(temp)
        #     temp = []

        # calc = readoutfilesobj.atoms_bands_objects[0].calc
        calc = readoutfilesobj.atoms_objects[0].calc
        # print(calc)
        dos = DOS(calc, width=0.2)
        d = dos.get_dos()
        e = dos.get_energies()
        # print(d)
        # print(e)


        # This is from the ASE website as an example
        # plt.plot(e, d)
        # plt.xlabel('energy [eV]')
        # plt.ylabel('DOS')
        # plt.show()

        self.y_to_plot.append(d)
        self.x_to_plot.append(e)
        self.labels.append(label)

if __name__ == "__main__":
    # You dont really need this. For testing we have left this block here.
    myplot = QEBandsReader("Sn.bands.out")
    path = kPathCreator()
    path.add_startval(G)
    path.add_kpath(X, 20)
    path.add_kpath(W, 20)
    path.add_kpath(L, 20)
    path.add_kpath(G, 20)
    path.add_kpath(K, 20)
    path.add_kpath(L, 20)
    path.out_kpath()
    # myplot.plot(path.k_path, path.highSymPoints_symbols)
    # myplot.calculate_mass(14,79,81)

    plotter= BandPlotter(path.k_distance, myplot.data, path.highSymPoints_position, path.highSymPoints_symbols)
    plotter.title = "Sn QE Band Diagram"
    # plotter.file_name = "Ef9shift"
    plotter.vlines = True
    plotter.set_y_range = False
    plotter.same_band_colour = True
    plotter.Ef_shift = 9
    plotter.plot()

    plotter2= BandPlotter(path.k_distance, myplot.data, path.highSymPoints_position, path.highSymPoints_symbols)
    plotter2.title = "Sn QE Band Diagram"
    # plotter2.file_name = "Ef9shift"
    plotter2.vlines = True
    plotter2.set_y_range = False
    plotter2.same_band_colour = True
    plotter2.Ef_shift = 0
    plotter2.band_colour = "r"
    # plotter2.new_fig = True
    plotter2.plot()