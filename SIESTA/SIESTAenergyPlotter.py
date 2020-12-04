"""
This file is only used to plot energy levels and not band diagrams.
The idea is to quickly read in sieta .EIG files for Gamma point calculations and plot them.
Multiple diagrams can be imported and then plotted together.
    This means importing multiple .EIG files and then pinning the diagrams to some level usually the vacuum level.
    This pinning level can be set manually or will be automatically read from the .out file.
"""

from ebk.SIESTA.SIESTAOutFileReader import SiestaReadOut
import matplotlib.pyplot as plt
import os
# print(os.getcwd())

class PlotEnergy():
    def __init__(self, *args, **kwargs):
        """
        The first argument is the file
        """
        self.file_name = "Energy_level_diagrams"
        self.SystemLabels = []
        self.Energies = []
        self.pin_levels = []
        self.line_widths = 0.2
        self.number_of_bands_to_plot = 30
        self.show_fermi = False  # Setting this to true will show fermi level
        self.show_homo = True
        self.show_lumo = True
        self.kwargs_list = []
        self.fermi_energies = []
        self.HOMOs = []
        self.LUMOs = []
        self.hlines = kwargs.get("hlines", True)
        self.vlines = kwargs.get("vlines", True)
        self.title = "Energy Levels"
        self.set_y_range_upper = kwargs.get("set_y_range_upper", True)
        self.set_y_range_lower = kwargs.get("set_y_range_lower", False)
        self.set_x_range = False
        self.Ef_shift = True
        self.y_to_plot = []
        self.x_to_plot = []
        self.labels = []
        self.same_band_colour = kwargs.get("same_band_colour", True)
        self.band_colour = ["b", "g", "r", "c", "m", "y", "k"]
        self.new_fig = False
        self.k_locations = None
        self.k_symbols = None
        self.dots = kwargs.get("dots", False)
        self.ylim_low = -9
        self.ylim_high = 0

    def load(self, file_name, pin_level, label = None, *args, **kwargs):
        """
        Will load the relevent files and values and prep for plotting
        inputs:
            file_name: Path to .EIG file
            pin_level: Level to which we should pin the energy levels
        """
        self.labels.append(label)
        self.pin_levels.append(pin_level)

        # First we will figure out the SystemLabel so we can see what we will be there
        SystemLabel = file_name.strip().split("/")[-1].split(".")[-2]
        path_to_file_components = file_name.strip().split("/")
        path_to_file = path_to_file_components[0]
        for x in range(1, len(path_to_file_components)-1):
            path_to_file = f"{path_to_file}/{path_to_file_components[x]}"

        # Failed attempt to get vac level automatically
        # print(f"./{path_to_file}/{SystemLabel}")
        # print(os.getcwd())
        # vac = SiestaReadOut(f"./{path_to_file}/{SystemLabel}").read_vacuum()
        # print(vac)
        # print(vac.file)
        # vac_max = vac.Vac_max

        # print(SystemLabel)
        # self.SystemLabels.append(SystemLabel)

        print(f"Opening file: {file_name}")
        file_lines = []
        with open(file_name, "r+") as file:
            for line in file:
                # print(line)
                file_lines.append(line)

        # Getting Fermi energies
        Ef=float(file_lines[0])
        self.fermi_energies.append(Ef-pin_level)

        temp_energies = []
        # This is the special case of the first data line that has a "1" in front.
        add_this = file_lines[2].split()
        for x in range(1,len(add_this)):
            temp_energies.append(float(add_this[x]))

        for x in range(3,len(file_lines)):
            add_this = file_lines[x].split()
            for i in range(len(add_this)):
                temp_energies.append(float(add_this[i]))
        temp_energies = [x - pin_level for x in temp_energies]
        # print(temp_energies)
        self.Energies.append(temp_energies)

        # Finding HOMO and LUMO levels
        old=temp_energies[0]
        new=temp_energies[0]
        for x in temp_energies:
            old=new
            new = x
            if new >= Ef:
                break
        self.HOMOs.append(old)
        self.LUMOs.append(new)
        # print(Ef, old, new)
        self.kwargs_list.append(kwargs)


    def plot(self, *args, **kwargs):
        """
        Plots all energy levels
        """
        fig, ax1 = plt.subplots(figsize=(4*(len(self.Energies)*self.line_widths*2+0.5), 10))
        # plt.rcParams["figure.figsize"] = (100,4)
        # plt.figure()
        # Setting y ranges
        if self.set_y_range_upper == True:
            ax1.set(ylim=(self.ylim_low,self.ylim_high))
            # ax1.set_ylim(self.ylim_high,None)
        # if self.set_y_range_lower == True:
        #     ax1.set_ylim(bottom=self.ylim_low)

        E = range(1,len(self.Energies)+1)
        for x in E:
            print(f"Plotting {x} of {len(E)}")
            for y in self.Energies[x-1]:
                # print(f".", end = "")
                x_coordinates = [x-self.line_widths, x+self.line_widths]
                y_coordinates = [y, y]
                plt.plot(x_coordinates, y_coordinates, **self.kwargs_list[x-1])
            if self.show_fermi:
                plt.plot(x_coordinates, [self.fermi_energies[x-1], self.fermi_energies[x-1]], "b--")
            if self.show_homo:
                plt.plot([x_coordinates[1]+0.1, x_coordinates[1]+0.1], [self.HOMOs[x-1], self.HOMOs[x-1]], "<b")
            if self.show_lumo:
                plt.plot([x_coordinates[1]+0.1, x_coordinates[1]+0.1], [self.LUMOs[x-1], self.LUMOs[x-1]], "<g")
        ax1.set_xticks(E)
        ax1.set_xticklabels(self.labels)
        ax1.set_xlabel("Ligand")
        ax1.set_ylabel("Energy (eV) (vac = 0)")
        ax1.set_title(f"{self.title}")
        plt.savefig(f"{self.file_name}.pdf")
        plt.show()



