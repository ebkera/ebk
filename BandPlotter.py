import matplotlib.pyplot as plt
import matplotlib
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


