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
        self.occupied_shifts = []
        self.unoccupied_shifts = []
        self.EAs = []
        self.IEs = []
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
        if "x_ticks_rotation" in kwargs:
            self.x_ticks_rotation = kwargs.get("x_ticks_rotation", 0)

    def load(self, file_name, pin_level, label = None, *args, **kwargs):
        """
        Will load the relevent files and values and prep for plotting
        inputs:
            file_name: Path to .EIG file (Should not contatain '.'s except for the one for thefile extension)
            pin_level: Level to which we should pin the energy levels. This will have to be input everytime.
        """
        self.labels.append(label)
        self.pin_levels.append(pin_level)

        # First we will figure out the SystemLabel so we can see what we will be there
        SystemLabel = file_name.strip().split("/")[-1].split(".")[-2]
        path_to_file_components = file_name.strip().split("/")
        path_to_file = path_to_file_components[0]
        for x in range(1, len(path_to_file_components)-1):
            path_to_file = f"{path_to_file}/{path_to_file_components[x]}"

        # Abandoned attempt at trying to do all steps for teh adjustment of energy levels from DSCF here.
        # path_to_file_NP1 = kwargs.get({"path_to_file_NP1":None})
        # path_to_file_NM1 = kwargs.get({"path_to_file_NM1":None})

        # if path_to_file_NM1 != None and path_to_file_NP1 != None:
        #     N = SiestaReadOut()
        #     NM1 = SiestaReadOut(path_to_file_NM1)
        #     NP1 = SiestaReadOut(path_to_file_NP1)
        #     self.E_NM1 = NM1.read_total_energy()
        #     self.E_NP1 = NP1.read_total_energy()
        #     self.E_N = N.read_total_energy()

        #     # print(ligand)
        #     # print(f"{(E_N)} {(E_NM1)} {(E_NP1)}")
        #     # if type(E_N) == None: NC = 0
        #     # if type(E_NM1) == None: E_NM1 = 0
        #     # if type(E_NP1) == None: print("haha")

        #     EA = -(E_NP1 - E_N)
        #     IE = (E_NM1 - E_N)
        #     Gap = E_NP1 + E_NM1 - 2*E_N

        #     plt.plot(NP_list.index(NP)+1, Gap, "bx", label = "Gap")
        #     leg_EA = plt.plot(NP_list.index(NP)+1, EA, "g<", label = "EA")
        #     leg_IE = plt.plot(NP_list.index(NP)+1, IE, "r<", label = "IE")

            # Failed attempt to get vac level automatically
            # print(f"./{path_to_file}/{SystemLabel}")
            # print(os.getcwd())
            # vac = SiestaReadOut(f"./{path_to_file}/{SystemLabel}").read_vacuum()
            # print(vac)
            # print(vac.file)
            # vac_max = vac.Vac_max

            # print(SystemLabel)
            # self.SystemLabels.append(SystemLabel)

        # Here we are planning to adjust the energy levels for the EA and IE that we got from DeltaSCF
        self.EAs.append(kwargs.get("EA",None))
        self.IEs.append(kwargs.get("IE",None))
        if "EA" in kwargs: del kwargs["EA"]
        if "IE" in kwargs: del kwargs["IE"]

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

        # Finding HOMO and LUMO levels
        old=temp_energies[0]
        new=temp_energies[0]
        for x in temp_energies:
            old=new
            new = x
            if new >= Ef - pin_level:
                break
        self.HOMOs.append(old)
        self.LUMOs.append(new)
        # print(Ef, old, new)
        self.kwargs_list.append(kwargs)

        # Adjusting the HOMO and LUMO levels if IE and EA are present:
        if self.EAs[-1] != None and self.IEs[-1] !=None:
            print(f"Adjusting occ/unocc levels as per DeltaSCF calculations")
            occupied_shift = (-1*self.IEs[-1]) - old  # Converting the postive IE to negative. adn old corresponds to homo
            unoccupied_shift = new - (-1*self.EAs[-1])
            # print(f"Delta_occ:\t{occupied_shift}") # For debugging
            # print(f"Delta_uoc:\t{unoccupied_shift}") # For debugging
            temp_energies = [x+occupied_shift if x<= Ef-pin_level else x-unoccupied_shift if x>=Ef-pin_level else x for x in temp_energies]
            # temp_energies = [x+occupied_shift for x in temp_energies if x<= Ef-pin_level]
            # temp_energies = [x-unoccupied_shift for x in temp_energies if x>= Ef-pin_level]
            self.HOMOs[-1]=self.HOMOs[-1]+occupied_shift
            self.LUMOs[-1]=self.LUMOs[-1]-unoccupied_shift
            # self.fermi_energies = self.fermi_energies[-1]

        # We are doing this at the end so we can adjust all the adjustments
        self.Energies.append(temp_energies)


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
            print(f"Plotting {x}\tof {len(E)}")
            for y in self.Energies[x-1]:
                # print(f".", end = "")
                x_coordinates = [x-self.line_widths, x+self.line_widths]
                y_coordinates = [y, y]
                plt.plot(x_coordinates, y_coordinates, **self.kwargs_list[x-1])
            if self.show_fermi:
                plt.plot(x_coordinates, [self.fermi_energies[x-1], self.fermi_energies[x-1]], "b--", label="Fermi Level")
            if self.show_homo:
                plt.plot([x_coordinates[1]+0.1, x_coordinates[1]+0.1], [self.HOMOs[x-1], self.HOMOs[x-1]], "<b", label="HOMO")
            if self.show_lumo:
                plt.plot([x_coordinates[1]+0.1, x_coordinates[1]+0.1], [self.LUMOs[x-1], self.LUMOs[x-1]], "<g", label="LUMO")
        ax1.set_xticks(E)
        ax1.set_xticklabels(self.labels)
        ax1.set_xlabel("Ligand")
        ax1.set_ylabel("Energy (eV) (vac = 0)")
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys(), loc = 'upper right')
        if hasattr(self, "x_ticks_rotation"): plt.xticks(rotation=self.x_ticks_rotation, ha='right')
        ax1.set_title(f"{self.title}")
        plt.savefig(f"{self.file_name}.pdf")
        plt.show()



