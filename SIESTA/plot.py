import os
import matplotlib
import sys
import subprocess
import matplotlib.pyplot as plt
import numpy as np
import shutil

# getting file paths
import pathlib
fpath_gnubands = pathlib.Path(__file__).parent.absolute()
print(fpath_gnubands)

def latexify(inputlist):
    """This method will convert a list of strings into it's latex form"""
    for index, item in enumerate(inputlist):
        if item == "Gamma":
            inputlist[index] = '$\Gamma$'
    return inputlist

class Band():
    def __init__(self, name):
        self.name = name
        self.file_name = name + ".bands.gnuplot.dat"
        self.bands_file_name = name + ".bands"
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

    def use_gnuband(self, command = "f"):
        """
        Executes the gnuband script and makes the SYSTEMLABEL.gnuband.dat file required for band plotting
        command: string (set of chars as commands)
            f or F: shift to fermi level
        """
        command = command.upper()
        if "F" in command:
            os.system(f"wsl {fpath_gnubands}\gnubands -F <  {self.file_name}.bands > {self.file_name}.bands.gnuplot.dat")
            print("Plot: Executing gnuplot for adjustment of values: Done")
        print(f"{fpath_gnubands}\gnubands")
        # os.system(f"{fpath_gnubands}\gnubands  < {SystemLabel}.bands > {SystemLabel}.bands.gnuplot.dat")

    def open_file(self):
        """Opens the file and gets all the data and also values. Full list of extractions include,
        From the gnuplot file - data and fermi energy
        From the .bands file - K path"""
        SystemLabel = self.name
        # print(fpath_gnubands)
        # os.system(f"{fpath_gnubands}\gnubands  < {SystemLabel}.bands > {SystemLabel}.bands.gnuplot.dat")

        file = open(self.file_name, "r")
        rows = []
        hashlinecount = 0  # looking at commented lines in order to extract fermi energy
        for line in file:
            rows.append(line)

        # Getting information from the comments section
        xdatatemp = []
        ydatatemp = []
        for row in rows:
            if "#" in row:  # Getting information from the comments
                hashlinecount += 1
                if hashlinecount == 7:  # Getting the E_f value
                    temp = row.split()
                    if len(temp) == 4:
                        self.E_f = temp[3]
                    elif len(temp) == 7:
                        self.E_f = temp[6]
                        self.Efshift = True
                    continue
                elif hashlinecount == 9:  # Getting the min and max band energies
                    temp = row.split()
                    self.E_min = float(temp[4])
                    self.E_max = float(temp[5])
                    continue
                elif hashlinecount == 11:  # Gettiing the minimum and maximum band numbers to plot
                    temp = row.split()
                    self.band_number_min = int(temp[5])
                    self.band_number_max = int(temp[6])
                    continue

            elif row.strip() == "":  # saving individual bands the first save will be a null band
                self.band_data_x.append(xdatatemp)  # Saving in a 2d list.
                self.band_data_y.append(ydatatemp)
                xdatatemp = []
                ydatatemp = []
            else:
                t = row.split()
                xdatatemp.append(float(t[0]))
                ydatatemp.append(float(t[1]))

        #  Getting rid of empty lists
        self.band_data_x = [x for x in self.band_data_x if x]
        self.band_data_y = [x for x in self.band_data_y if x]

        # shifting lines by E_f
        # self.band_data_y_Efshift = self.band_data_y + self.E_f

        #  opening the bands file to get the kpath
        bands_file = open(self.bands_file_name, "r")
        bands_file_rows = []
        hashlinecount = 0  # looking at commented lines in order to extract fermi energy
        for bands_file_line in bands_file:
            bands_file_rows.append(bands_file_line)

        length_bands_file_rows = len(bands_file_rows)
        still_two = True  # Check to see if there are only two things in the line
        for i in range(length_bands_file_rows-1, 0, -1):
            temperary = bands_file_rows[i].split()
            if len(temperary) == 2:
                self.kpath_k.append(float(temperary[0].strip()))
                self.kpath_symbol.append(str(temperary[1].strip("'")))
            elif len(temperary) != 2:
                break
        self.kpath_symbol = latexify(self.kpath_symbol)

    def plotter(self, *args, **kwargs):
        bands_to_plot = 0
        if len(args) != 0:
            bands_to_plot = []
            bands_to_plot_temp = args[0]  # First argument
            for x in bands_to_plot_temp:
                bands_to_plot.append(int(x))

        length = len(self.band_data_x)  # Number of original bands, have to select from this to plot
        plt.figure()
        # Plotting vertical lines so that we know the exact positions of the high symmetry points
        if kwargs["h"]:
            for xc in self.kpath_k:  # Plotting a dotted line that will at high symmetry points on the Kpath
                plt.axvline(x=xc, linestyle='-', color='k', linewidth=0.1)

        # Plotting a horizontal line so that we know the Fermi Energy level
        plt.axhline(linewidth=0.1, color='k')

        # The plotting stuff
        if bands_to_plot == 0:
            for i in range(1, length):
                x = self.band_data_x[i]
                y = self.band_data_y[i]
                plt.plot(x, y, linewidth=0.4)
        else:
            for band in bands_to_plot:
                x = self.band_data_x[band]
                y = self.band_data_y[band]
                plt.plot(x, y, linewidth=0.4)

        plt.title("Band diagram for bulk " + self.name + " (E$_f$ = " + str(self.E_f) +"eV)")
        # plt.plot([self.band_data_x[1][0], self.band_data_x[1][len(self.band_data_x[1]) - 1]], [self.E_f, self.E_f],
        #          '--k', label="E$_f$")  # Plotting the E_f
        plt.xlabel('k ($\\frac{\\pi}{a}$)')
        plt.ylabel('Energy in eV (E - E$_f$)')
        # plt.legend(loc='upper left')
        # plt.yticks(np.linspace(-23,150,11) , np.linspace(-23,150,11))
        # plt.text(0.8, 0.8, r'E$_f$ = ',fontsize=20)
        plt.xticks(self.kpath_k, self.kpath_symbol)


        # trying to write to file
        write_test = False  # This becomes True if the image is saved
        counter = 0  # This keeps track of the No. of attempts to save the figure
        try:
            plt.savefig(self.name + "_BandDiagramPlots.pdf")
            print("*** Saved figure with default file name")
        except:
            while write_test is False and counter < 5:
                try:
                    plt.savefig(self.name + "_BandDiagramPlots" + str(counter) + ".pdf")
                    print("!** Could not save figure with default file name. "
                          "Probably file is open and permissions are denied")
                    print("*** Saved figure as numbered version of default name")
                    write_test = True
                except:
                    counter += 1
            if write_test is False and counter == 5:
                print("!** Too many attempts to save as numbered version of default file name."
                      " Suggested to close program that is denying write access.")

    def execute_gnuplot(self, shift_fermi = True):
        # os.system("./gnubands -F < " + SystemLabel+".bands > " + SystemLabel + ".bands.gnuplot.dat")
        print(os.getcwd())
        path = self.name.split("/")
        print(path)
        newpath = "."
        for x in range(1, len(path)-1):
            newpath = f"{newpath}/{path[x]}"
        print(f"new: {newpath}")
        shutil.copy(f"{fpath_gnubands}\gnubands", f"{newpath}/gnubands")
        print(f"{newpath}/gnubands")
        print(f"{self.name}.bands.gnuplot.dat")
        print(f"{self.name}.bands")
        os.system(f"wsl {newpath}/gnubands -F < {self.name}.bands > {self.name}.bands.gnuplot.dat")
        print("*** Finished executing gnuplot for adjustment of values")

if __name__ == "__main__":
    SystemLabel = "Sn"
    arg_list_length = len(sys.argv)
    bands_to_plot_main = [x for x in sys.argv if x.isdigit()]  # getting the bands to plot
    if "r" in sys.argv:  # r for "run seista"
        os.system("siesta < " + SystemLabel+".fdf | tee Sn.out ")
    if "g" in sys.argv:  # g for using gnubands file this will shift the values
        print("*** Starting to execute gnuplot for adjustment of values")
        os.system("./gnubands -F < " + SystemLabel+".bands > " + SystemLabel + ".bands.gnuplot.dat")
        print("*** Finished executing gnuplot for adjustment of values")
    if "d" in sys.argv:  # d for "don't shift values when using gnubands file"
        print("*** Starting to execute gnuplot for adjustment of values - No shifting done")
        os.system("./gnubands  < " + SystemLabel + ".bands > " + SystemLabel + ".bands.gnuplot.dat")
        print("*** Finished executing gnuplot for adjustment of values")

    # bands.plotter(6, 9, 10, 11, 12)
    bands = Band(SystemLabel)
    bands.open_file()

    if len(bands_to_plot_main) != 0:
        print("***Plotting custom band set:")
        print(bands_to_plot_main)
        if 'h' in sys.argv:
            bands.plotter(bands_to_plot_main, h=True)
        else:
            bands.plotter(bands_to_plot_main, h=False)

    else:
        if 'h' in sys.argv:
            bands.plotter(h=True)
        else:
            bands.plotter(h=False)
