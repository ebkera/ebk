import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from ase.dft.kpoints import resolve_custom_points, find_bandpath_kinks
from ase.dft.dos import DOS
from matplotlib import gridspec

class BandPlotter_old():
    """Retained as legacy code Please use the below code"""
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
        # self.number_of_bands_to_plot = 30
        self.fermi_level = 0
        self.hlines = False
        self.vlines = False
        self.title = "Band diagram"
        self.set_y_range = False
        self.ylim_low = -5
        self.ylim_high = 5
        self.xlim_low = 0
        self.xlim_high = 0
        self.Ef_shift = True  # notimplemented yet
        self.y_to_plot = []
        self.same_band_colour = False
        self.band_colour = "b"
        self.new_fig = False
        self.x_margins = 0
        self.saveas_extension = "pdf"
        self.extra_horizontal_lines = []  # [ycoord,label]
        self.extra_vertical_lines = []
        self.plot_arrow = [0,0,0,0]
        self.include_xlabel = True
        self.include_ylabel = True
        self.include_title = True

    def plot(self):
        """
        |All the features of the band plot are set here
        |Inputs: None
        """
        for i in range(0,len(self.y)):
            self.y_to_plot.append([x - self.fermi_level for x in self.y[i]])
        
        # To get multiple band diagrams together just plot them on the same figure (with the same file name) and when ever you want a new figure with just the new plots
        #  just do new_fig = True. If you want to save individual figures just give a new file name (with new_fig = True).
        if self.new_fig == True:
            plt.figure()
        
        # Setting vertical lines
        if self.vlines == True:
            for xc in self.k_locations:  # Plotting a dotted line that will run vertical at high symmetry points on the Kpath
                plt.axvline(x=xc, linestyle='dotted', color='k', linewidth=0.05)

        # Setting horizontal lines Usualy only fermi
        if self.hlines == True:
            plt.axhline( linestyle='dotted',linewidth=0.05, color='k')

        # Setting y ranges
        if self.set_y_range == True:
            plt.ylim([self.ylim_low,self.ylim_high])

        # We plot the figure here
        for band in range(0,len(self.y)):
            if self.same_band_colour == True:
                plt.plot(self.x, self.y_to_plot[band], self.band_colour)
            else:
                plt.plot(self.x, self.y_to_plot[band])

        if len(self.extra_horizontal_lines) !=0:
            for line in self.extra_horizontal_lines:
                # print(line)
                plt.axhline(y=line[0], label = f"{line[1]}")

        plt.xticks(self.k_locations, self.k_symbols)
        if self.include_xlabel: plt.xlabel("$\\vec{{k}}$")
        if self.include_ylabel: plt.ylabel("E-E$_f$ (eV)")
        if self.include_title: plt.title(f"{self.title}")
        if self.plot_arrow != [0,0,0,0]:
            # kwargs = {}
            # kwargs.update({"width":0.00})
            # kwargs.update({"linestyle":"--"})
            # kwargs.update({"length_includes_head":True})
            # kwargs.update({"head_width":0.0001})
            # kwargs.update({"head_length":6})
            # plt.arrow(self.plot_arrow[0],self.plot_arrow[1],self.plot_arrow[2],self.plot_arrow[3], **kwargs)
            plt.annotate("", xy=(self.plot_arrow[2], self.plot_arrow[3]), xytext=(self.plot_arrow[0], self.plot_arrow[1]), arrowprops=dict(arrowstyle="->"))
        # plt.arrow(2,2,2,2)
        # width1 = 5
        # height1 = 5
        # width_height_1 = (width1, height1)
        # plt.figure(figsize=width_height_1)
        plt.margins(x=self.x_margins)
        # plt.tight_layout()
        # plt.legend(loc= "upper right")
        plt.savefig(f"Bands_{self.file_name}.{self.saveas_extension}")
        plt.show()

class BandPlotter():
    def __init__(self, **kwargs):
        """
        If plotting multiple band plots here can do only same path.
        The idea of how to use this class is as below:
            1) you set the required parameters that you need to plot the graph
            2) You add all the bands you want to add by using the add_to_plot() method
            3) Plot all bands that you have added using above method by calling the plot() method
        """
        self.file_name = kwargs.get("file_name", "Band diagram")
        self.title = kwargs.get("title", self.file_name)
        self.number_of_bands_to_plot = 30
        self.fermi_level = 0
        self.hlines = kwargs.get("hlines", True)
        self.vlines = kwargs.get("vlines", True)
        self.dos_title = "Density of States"
        self.set_y_range = kwargs.get("set_y_range", True)
        self.set_x_range = False
        self.Ef_shift = True
        self.y_to_plot = []
        self.x_to_plot = []
        self.dos = []
        self.E_dos = []
        self.dos_color_fills = []
        self.labels = []
        self.same_band_colour = kwargs.get("same_band_colour", True)
        self.band_colour = ["b", "g", "r", "c", "m", "y", "k"]
        self.band_colour = ["blue","red", "green", "purple", "brown", "pink", "orange", "gray", "olive", "cyan", "blue", "red", "green", "purple", "brown", "pink", "orange", "gray", "olive", "cyan"]
        self.new_fig = False
        self.k_locations = None
        self.k_symbols = None
        self.dots = kwargs.get("dots", False)
        self.include_dos = kwargs.get("include_dos", False)
        self.include_tdm = kwargs.get("include_tdm", False)  # include transitional dipole matrix
        self.dos_units = kwargs.get("dos_units", "states/eV")
        self.plot_only_dos = kwargs.get("only_dos", False)
        self.pin_fermi = kwargs.get("pin_fermi", "scf") # if there are difference in fermi levels (Eg E_scf and E_nscf) then use this to pin the levels There options "off", "scf", "nscf" 
        self.plot_from_QE_dos_file = True  # This will make the code read in files from a QE dos output file.
        self.QE_dos_file_path = ""
        self.ylim_low = -2
        self.ylim_high = 2
        self.x_margins = 0
        self.xlim_low = 0
        self.xlim_high = 300
        self.plt_width = 3.5
        self.plt_height = 3.5
        self.saveas_extension = "pdf"
        self.include_bandgaparrow = True
        self.arrow_data = kwargs.get("arrow_data", [])
        self.show_figs = True
        self.legend = True
        self.dos_ylabel = ""
        self.extra_data_x = []
        self.extra_data_y = []
        self.extra_kwargs = []
        self.extra_labels = []
        self.tdm_k        = []
        self.tdm_tdm      = []
        self.tdm_labels   = []

    def plot(self):
        """
        |All the features of the band plot are set here
        |Inputs: None
        """
        # Setting the dimensions of the saved image
        if self.include_dos:
            fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, constrained_layout=True, figsize=(self.plt_width,self.plt_height))
            gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1], wspace=0)
            ax1 = plt.subplot(gs[0])
            ax2 = plt.subplot(gs[1])
            ax2.yaxis.tick_right()
            ax2.yaxis.set_label_position("right")
            ax2.xaxis.get_major_ticks()[0].draw = lambda *args:None
        elif self.include_tdm:
            fig, (ax1, ax2) = plt.subplots(2, 1, constrained_layout=True, figsize=(self.plt_width,self.plt_height))
            gs = gridspec.GridSpec(2, 1,height_ratios=[3, 1], hspace=0.3)
            ax1 = plt.subplot(gs[0])
            ax3 = plt.subplot(gs[1])
        elif self.plot_only_dos:
            fig, ax2 = plt.subplots(figsize=(self.plt_width,self.plt_height))
        else:
            fig, ax1 = plt.subplots(figsize=(self.plt_width,self.plt_height))

        # #here we do an alternative
        # fig, ax = plt.subplots(2, 2, sharey=True, constrained_layout=True, figsize=(self.plt_width,self.plt_height))
        # if self.include_dos:
        #     gs = gridspec.GridSpec(2, 2, width_ratios=[3, 1], wspace=0, hspace=0)
        #     ax1 = plt.subplot(gs[0][0])
        #     ax2 = plt.subplot(gs[0][1])
        #     ax3 = plt.subplot(gs[1][0])
        #     ax2.yaxis.tick_right()
        #     ax2.yaxis.set_label_position("right")
        #     ax2.xaxis.get_major_ticks()[0].draw = lambda *args:None
        # elif self.plot_only_dos:
        #     fig, ax2 = plt.subplots(figsize=(self.plt_width,self.plt_height))
        # else:
        #     fig, ax1 = plt.subplots(figsize=(self.plt_width,self.plt_height))


        # Setting vertical lines
        if self.vlines == True:
            if not self.plot_only_dos:
                for xc in self.k_locations:  # Plotting a dotted line that will run vertical at high symmetry points on the Kpath
                    ax1.axvline(x=xc, linestyle='--', color='k', alpha=0.2)
                    # try: ax3.axvline(x=xc, linestyle='--', color='k', alpha=0.2)
                    # except: pass
        # Setting horizontal line on the fermi energy
        if self.hlines == True:
            if self.plot_only_dos or self.include_dos:
                ax2.axhline(linewidth=0.1, color='k', alpha=0.2)
            if not self.plot_only_dos :
                ax1.axhline(linewidth=0.1, color='k', alpha=0.2)
                # try: ax3.axhline(linewidth=0.1, color='k', alpha=0.2)
                # except: pass
        # Setting y ranges
        if self.set_y_range == True:
            if self.plot_only_dos:
                ax2.set_ylim([self.ylim_low, self.ylim_high])
            elif self.include_dos:
                ax1.set_ylim([self.ylim_low, self.ylim_high])
                ax2.set_ylim([self.ylim_low, self.ylim_high])
            else:
                ax1.set_ylim([self.ylim_low, self.ylim_high])
        if self.set_x_range == True:
            if self.plot_only_dos or self.include_dos:
                # This for now applies only to DOS
                ax2.set_xlim([self.xlim_low, self.xlim_high])
            if not self.plot_only_dos :
                pass
        # we plot the dos figure here
        if self.plot_only_dos or self.include_dos:
            for i,v in enumerate(self.dos):
                list_of_zeros = [0] * len(self.dos[i])
                ax2.plot(self.dos[i], self.E_dos[i], self.band_colour[i], label = self.labels[i])
                if self.dos_color_fills[i] == True: 
                    ax2.fill_betweenx(self.E_dos[i], 0, self.dos[i], facecolor=self.band_colour[i], alpha = 0.2)
                    # ax2.fill_betweenx(self.dos[i], self.E_dos[i], facecolor=self.band_colour[i], alpha = 0.2)

            ax2.set_xlabel(self.dos_units)
            ax2.set_ylabel(self.dos_ylabel)
            ax2.set_title(f"{self.dos_title}")
            ax2.tick_params(axis="y",direction="in")
            ax2.tick_params(axis="x",direction="in")
            ax2.axes.yaxis.set_ticklabels([])

        if self.include_tdm:
            for i in range(len(self.tdm_k)):
                ax3.fill_between(self.tdm_k[i], self.tdm_tdm[i], label=self.tdm_labels[i], alpha = 0.2)
            # ax3.set_xticks(self.k_locations)
            # ax3.set_xticklabels(self.k_symbols)
            ax3.set_xticklabels([])
            ax3.set_xticks([])
            ax3.set_ylabel(f"P$^2$ (D$^2$)")
            ax3.set_title(f"")
            ax3.tick_params(axis="y",direction="in")
            # ax3.tick_params(axis="x",direction="in")
            # ax3.set_xlabel("$\\vec{{k}}$")
            ax3.margins(x=self.x_margins)


            # ax2.legend(loc="upper right")
        if not self.plot_only_dos :
            # We plot the bands figure here
            for i,v in enumerate(self.y_to_plot):
                for band in v:
                    if self.same_band_colour == True:
                        if self.labels[i] != None:
                            ax1.plot(self.x_to_plot[i], band, self.band_colour[i], label = self.labels[i],linewidth=0.5)
                            self.labels[i] = None  # Here we have set the label to none since we are plotting multiple bands with teh same label and we dont want to list multiple labels
                        else:
                            ax1.plot(self.x_to_plot[i], band, self.band_colour[i], label = self.labels[i],linewidth=0.5)
                    else:
                        ax1.plot(self.x_to_plot[i], band)  # Here we have ignored labels since individual bands have no labels

            # if not self.include_tdm:
            ax1.set_xticks(self.k_locations)
            ax1.set_xticklabels(self.k_symbols)
            # else: 
            #     ax1.set_xticks([])
            #     ax1.set_xticklabels([])
            ax1.tick_params(axis="y",direction="in")
            ax1.tick_params(axis="x",direction="in")
            ax1.set_xlabel("$\\vec{{k}}$")
            ax1.set_ylabel("Energy (eV)")
            ax1.set_title(f"{self.title}")
            ax1.margins(x=self.x_margins)
            if self.include_bandgaparrow and self.arrow_data != []:
                for i,arrow in enumerate(self.arrow_data):
                    ax1.annotate('', xytext=(arrow[0],arrow[1]), xy=(arrow[2],arrow[3]), arrowprops={'arrowstyle': '->', 'color':self.band_colour[i]}, va='center', color='red')
            for i in range(len(self.extra_data_x)):
                ax1.plot(self.extra_data_x[i], self.extra_data_y[i], label=self.extra_labels[i], **self.extra_kwargs[i])
            if self.legend:  ax1.legend(loc="upper right")
        plt.margins(x=self.x_margins)
        plt.tight_layout()
        plt.savefig(f"Bands_{self.file_name}.{self.saveas_extension}")
        if self.show_figs: plt.show()

    def add_to_plot(self, Ef, k_dist, bands, hsp, hss, E_dos=0, dos=0, dos_shift=0, label=None, dos_color_fill=True):
        """
        |Here you add individual plots that need to be plot and then just plot them with the plot() method
        """

        if dos != 0:
            E_dos = [x-dos_shift for x in E_dos]
            self.dos.append(dos)
            self.E_dos.append(E_dos)
            self.dos_color_fills.append(dos_color_fill)

        # This is for the band structure
        self.k_locations = hsp
        self.k_symbols = hss
        Energy_to_plot = []
        if self.Ef_shift == True:
            for band in bands:
                if self.pin_fermi != "nscf":
                    # print(Ef)
                    Energy_to_plot.append([E - Ef for E in band])
                else:
                    Energy_to_plot.append([E - Ef_nscf for E in band])  # here the assumption is that if this option is ever reached then the idea is that an nscf calculation has already being done and dos is being plotted.
        self.y_to_plot.append(Energy_to_plot)
        self.x_to_plot.append(k_dist)
        self.labels.append(label)
    
    def add_to_extras(self, x, y, label, **kwargs):
        self.extra_data_x.append(x)
        self.extra_data_y.append(y)
        self.extra_labels.append(label)
        self.extra_kwargs.append(kwargs)

    def add_to_TDM(self, energies, tdm, tdm_labels=""):
        self.tdm_k.append(energies)
        self.tdm_tdm.append(tdm)
        self.tdm_labels.append(tdm_labels)

class BandPlotterASE():
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
        self.hlines = kwargs.get("hlines", True)
        self.vlines = kwargs.get("vlines", True)
        self.title = "Band diagram"
        self.dos_title = "Density of States"
        self.set_y_range = kwargs.get("set_y_range", True)
        self.set_x_range = False
        self.Ef_shift = True
        self.y_to_plot = []
        self.x_to_plot = []
        self.dos = []
        self.E_dos = []
        self.labels = []
        self.same_band_colour = kwargs.get("same_band_colour", True)
        self.band_colour = ["b", "g", "r", "c", "m", "y", "k"]
        self.new_fig = False
        self.k_locations = None
        self.k_symbols = None
        self.dots = kwargs.get("dots", False)
        self.include_dos = kwargs.get("include_dos", False)
        self.plot_only_dos = kwargs.get("only_dos", False)
        self.pin_fermi = kwargs.get("pin_fermi", "scf") # if there are difference in fermi levels (Eg E_scf and E_nscf) then use this to pin the levels There options "off", "scf", "nscf" 
        self.plot_from_QE_dos_file = True  # This will make the code read in files from a QE dos output file.
        self.QE_dos_file_path = ""
        self.ylim_low = -5
        self.ylim_high = 5
        self.x_margins = 0
        # self.xlim_low = 0
        # self.xlim_high = 0
        self.saveas_extension = "pdf"

    # def get_dos(self, readoutfilesobj):
    #     """
    #     This method gets the dos parameters for plotting
    #     This can be an ordinary function and does not have to be a method for this class but I have included it as such so that we can be more flexible in the future
    #     """
    #     calc = readoutfilesobj.atoms_nscf_objects[0].calc
    #     dos = DOS(calc, width=0.2)
    #     self.dos.append(dos.get_dos())
    #     self.E_dos.append(dos.get_energies())

    def plot(self):
        """
        |All the features of the band plot are set here
        |Inputs: None
        """
        # Setting the dimensions of the saved image
        plt.rcParams["figure.figsize"] = (8,4)
        if self.include_dos:
            # fig, (ax1, ax2) = plt.subplots(1,2)
            gs = gridspec.GridSpec(1, 2, width_ratios=[4, 1])
            gs.update(wspace=0.25, hspace=0.05)  #this does not seem to work.
            ax1 = plt.subplot(gs[0])
            ax2 = plt.subplot(gs[1])
        elif self.plot_only_dos:
            fig, ax2 = plt.subplots()
        else:
            fig, ax1 = plt.subplots()

        # Setting vertical lines
        if self.vlines == True:
            if not self.plot_only_dos:
                for xc in self.k_locations:  # Plotting a dotted line that will run vertical at high symmetry points on the Kpath
                    ax1.axvline(x=xc, linestyle='-', color='k', linewidth=0.1)
        # Setting horizontal line on the fermi energy
        if self.hlines == True:
            if self.plot_only_dos or self.include_dos:
                ax2.axhline(linewidth=0.1, color='k')
            if not self.plot_only_dos :
                ax1.axhline(linewidth=0.1, color='k')
        # Setting y ranges
        if self.set_y_range == True:
            if self.plot_only_dos:
                ax2.set_ylim([self.ylim_low, self.ylim_high])
            elif self.include_dos:
                ax1.set_ylim([self.ylim_low, self.ylim_high])
                ax2.set_ylim([self.ylim_low, self.ylim_high])
            else:
                ax1.set_ylim([self.ylim_low, self.ylim_high])
        if self.set_x_range == True:
            if self.plot_only_dos or self.include_dos:
                # This for now applies only to DOS
                ax2.set_xlim([self.xlim_low, self.xlim_high])
            if not self.plot_only_dos :
                pass
        # we plot the dos figure here
        if self.plot_only_dos or self.include_dos:
            # print(self.include_dos)
            for i,v in enumerate(self.y_to_plot):
                ax2.plot(self.dos[i], self.E_dos[i], self.band_colour[i], label = self.labels[i])
            ax2.set_xlabel("DOS")
            ax2.set_ylabel("Energy (eV)")
            ax2.set_title(f"{self.dos_title}")
            ax2.legend()
        if not self.plot_only_dos :
            # We plot the bands figure here
            for i,v in enumerate(self.y_to_plot):
                for band in v:
                    if self.same_band_colour == True:
                        if self.labels[i] != None:
                            ax1.plot(self.x_to_plot[i], band, self.band_colour[i], label = self.labels[i])
                            self.labels[i] = None  # Here we have set the label to none since we are plotting multiple bands with teh same label and we dont want to list multiple labels
                        else:
                            ax1.plot(self.x_to_plot[i], band, self.band_colour[i], label = self.labels[i])
                    else:
                        ax1.plot(self.x_to_plot[i], band)  # Here we have ignored labels since individual bands have no labels
            print(self.k_locations)
            print(self.k_symbols)
            ax1.set_xticks( self.k_locations)
            ax1.set_xticklabels(self.k_symbols)
            ax1.set_xlabel("K path")
            ax1.set_ylabel("Energy (eV)")
            ax1.set_title(f"{self.title}")
            ax1.legend()
        plt.margins(x=self.x_margins)
        plt.savefig(f"Bands_{self.file_name}.{self.saveas_extension}")
        plt.show()

    def add_to_plot(self, readoutfilesobj, label = None):
        """
        |Here you add individual plots that need to be plot and then just plot them with the plot() method
        |Use this method which is a part of the BandPlotterASE class you will have to give the bands to plot using a readoutfilesobj type of object
        """
        try:
            Ef = readoutfilesobj.atoms_objects[0].calc.get_fermi_level()
        except:
            print(f"add_to_plot: Could not read fermi level on scf file")
        try:
            kpts = readoutfilesobj.atoms_bands_objects[0].calc.get_ibz_k_points()
        except:
            print(f"add_to_plot: Could not read k points on bands file")

        if self.include_dos or self.plot_only_dos:
            if self.plot_from_QE_dos_file:
                E = []
                dos = []
                with open(readoutfilesobj.dos_file_name, "r") as dos_file:
                    for line in dos_file:
                        # print(line)
                        if "#" not in line:
                            data = line.strip().split()
                            E.append(float(data[0]) - Ef)
                            dos.append(float(data[1]))
                # self.y_to_plot.append(dos)
                # self.x_to_plot.append(E)
                self.dos.append(dos)
                self.E_dos.append(E)
                if self.plot_only_dos:  # Other wise we wil lappend to labels twice sicnce we will do it again at theloading of bands
                    self.labels.append(label)

            else:
                # we are here getting the required information for DOS
                # self.get_dos(readoutfilesobj)
                calc = readoutfilesobj.atoms_nscf_objects[0].calc
                Ef_nscf = readoutfilesobj.atoms_nscf_objects[0].calc.get_fermi_level()
                dos = DOS(calc, width=0.2)
                self.dos.append(dos.get_dos())
                E_nscf = dos.get_energies()
                Ef_diff = Ef - Ef_nscf
                if abs(Ef_diff) >= 0.001:
                    print(f"Warning!!! In {readoutfilesobj.identifier[0]}:\t|E_f_scf-E_f_nscf| = {abs(Ef_diff):.3} eV (>= 0.001 eV); E_f_scf = {Ef} eV,  E_f_nscf = {Ef_nscf} eV")
                if self.pin_fermi != "scf": self.E_dos.append(E_nscf)
                else: self.E_dos.append([E - Ef_diff for E in E_nscf])

        # Test space for k path and k high symmetry points
        # print(kpts)
        # print(self.plot_only_dos)
        if self.plot_only_dos == False:
            # This is for teh band structure
            path = readoutfilesobj.atoms_bands_objects[0].cell.bandpath(npoints=0)
            kinks = find_bandpath_kinks(readoutfilesobj.atoms_bands_objects[0].cell, kpts, eps=1e-5)  # These are the high symmetry points in use 
            pathspec = resolve_custom_points(kpts[kinks], path.special_points, eps=1e-5) # This gives the postions for the relevant high symmetry points
            # path.kpts = kpts
            # path.path = pathspec

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
                    if self.pin_fermi != "nscf":
                        Energy_to_plot.append([E - Ef for E in band])
                    else:
                        Energy_to_plot.append([E - Ef_nscf for E in band])  # here the assumption is that if this option is ever reached then the idea is that an nscf calculation has already being done and dos is being plotted.
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



# copy of the class before the implementation of all DOS related things into the things

# class BandPlotterASE():
#     # Under costruction!!
#     def __init__(self, **kwargs):
#         """
#         If plotting multiple band plots here can do only same path.
#         The idea of how to use this class is as below:
#             1) you set the required parameters that you need to plot the graph
#             2) You add all the bands you want to add by using the add_to_plot() method
#             3) Plot all bands that you have added using above method by calling the plot() method
#         """
#         self.file_name = "band_diagram"
#         self.number_of_bands_to_plot = 30
#         self.fermi_level = 0
#         self.hlines = False
#         self.vlines = False
#         self.title = "Band diagram"
#         self.set_y_range = False
#         self.ylim_low = -5
#         self.ylim_high = 5
#         self.xlim_low = 0
#         self.xlim_high = 0
#         self.Ef_shift = True
#         self.y_to_plot = []
#         self.x_to_plot = []
#         self.labels = []
#         self.same_band_colour = False
#         self.band_colour = ["b", "g", "r", "c", "m", "y", "k"]
#         self.new_fig = False
#         self.k_locations = None
#         self.k_symbols = None
#         self.dots = kwargs.get("dots", False)

#     def plot(self):
#         """
#         |All the features of the band plot are set here
#         |Inputs: None
#         """

#         # Setting the dimensions of the saved image
#         plt.rcParams["figure.figsize"] = (14,9)

#         # Setting vertical lines
#         if self.vlines == True:
#             for xc in self.k_locations:  # Plotting a dotted line that will run vertical at high symmetry points on the Kpath
#                 plt.axvline(x=xc, linestyle='-', color='k', linewidth=0.1)

#         # Setting horizontal line on the fermi energy
#         if self.hlines == True:
#             plt.axhline(linewidth=0.1, color='k')

#         # Setting y ranges
#         if self.set_y_range == True:
#             plt.ylim([self.ylim_low, self.ylim_high])

#         # We plot the figure here
#         for structure in range(0,len(self.y_to_plot)):
#             # print(f"structure: {}")
#             for band in range(0,len(self.y_to_plot[structure])):
#                 if self.same_band_colour == True:
#                     if self.labels[structure] != None:
#                         plt.plot(self.x_to_plot[structure], self.y_to_plot[structure][band], self.band_colour[structure], label = self.labels[structure])
#                         self.labels[structure] = None
#                     else:
#                         plt.plot(self.x_to_plot[structure], self.y_to_plot[structure][band], self.band_colour[structure], label = self.labels[structure])
#                 else:
#                     plt.plot(self.x_to_plot[structure], self.y_to_plot[structure][band])


#         # # This is only for testing with Pierre and should be removed once done
#         # # We plot the figure here
#         # for structure in range(0,len(self.y_to_plot)):
#         #     # print(f"structure: {}")
#         #     for band in range(0,len(self.y_to_plot[structure])):
#         #         if self.same_band_colour == True:
#         #             if self.labels[structure] != None:
#         #                 plt.plot(self.x_to_plot[structure], self.y_to_plot[structure][band], ".", label = self.labels[structure])
#         #                 self.labels[structure] = None
#         #             else:
#         #                 plt.plot(self.x_to_plot[structure], self.y_to_plot[structure][band], ".", label = self.labels[structure])
#         #         else:
#         #             plt.plot(self.x_to_plot[structure], self.y_to_plot[structure][band], ".")


#         plt.xticks(self.k_locations, self.k_symbols)
#         plt.xlabel("K path")
#         plt.ylabel("Energy (eV)")
#         plt.title(f"{self.title}")
#         plt.legend()
#         plt.savefig(f"Bands_{self.file_name}.pdf")
#         plt.show()

#     def add_to_plot(self, readoutfilesobj, label = None):
#         """
#         |Here you add individual plots that need to be plot and then just plot them with the plot() method
#         |Use this method which is a part of the BandPlotterASE class you will have to give the bands to plot using a readoutfilesobj type of object
#         """
#         try:
#             Ef = readoutfilesobj.atoms_objects[0].calc.get_fermi_level()
#         except:
#             print(f"add_to_plot: Could not read fermi level")
#         try:
#             kpts = readoutfilesobj.atoms_bands_objects[0].calc.get_ibz_k_points()
#         except:
#             print(f"add_to_plot: Could not read k points")

#         # Test space for k path and k high symmetry points
#         # print(kpts)
#         path = readoutfilesobj.atoms_bands_objects[0].cell.bandpath(npoints=0)
#         kinks = find_bandpath_kinks(readoutfilesobj.atoms_bands_objects[0].cell, kpts, eps=1e-5)  # These are the high symmetry points in use 
#         pathspec = resolve_custom_points(kpts[kinks], path.special_points, eps=1e-5) # This gives the postions for the relevant high symmetry points
#         path.kpts = kpts
#         path.path = pathspec

#         klengths = []
#         for x in range(0, len(kpts)):
#             if x == 0:
#                 kdist = np.sqrt((0.0 - kpts[x][0])**2 + (0.0 - kpts[x][1])**2 +(0.0 - kpts[x][2])**2)
#                 klengths.append(kdist)
#             else:
#                 kdist = np.sqrt((kpts[x-1][0] - kpts[x][0])**2 + (kpts[x-1][1] - kpts[x][1])**2 + (kpts[x-1][2]- kpts[x][2])**2)
#                 klengths.append(kdist+klengths[x-1])
#         self.k_locations = []

#         for x in range(len(kinks)):
#             self.k_locations.append(klengths[kinks[x]])

#         self.k_symbols = []
#         for x in pathspec:
#             if x == "G":
#                 self.k_symbols.append("$\Gamma$")
#             else:
#                 self.k_symbols.append(x)

#         energies = []
#         for s in range(readoutfilesobj.atoms_bands_objects[0].calc.get_number_of_spins()):
#             energies.append([readoutfilesobj.atoms_bands_objects[0].calc.get_eigenvalues(kpt=k, spin=s) for k in range(len(kpts))])
#         # print(f"lenght of kpoints: {range(len(kpts))}")
#         Energy_to_plot = []
#         if self.Ef_shift == True:
#             for band in energies[0]:
#                 Energy_to_plot.append([E - Ef for E in band])
#         tempMain = []
#         temp = []
#         for x in range(len(Energy_to_plot[0])):
#             for kpoint in Energy_to_plot:
#                 temp.append(kpoint[x])
#             tempMain.append(temp)
#             temp = []
#         self.y_to_plot.append(tempMain)
#         self.x_to_plot.append(klengths)
#         self.labels.append(label)
