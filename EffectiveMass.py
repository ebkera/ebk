import os
import matplotlib
matplotlib.use('Agg')  # no UI backend required if working in the wsl without a UI
import sys
import subprocess
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np
from matplotlib import gridspec

def latexify(inputlist):
    """This method will convert a list of strings into it's latex form"""
    for index, item in enumerate(inputlist):
        if item == "Gamma":
            inputlist[index] = '$\Gamma$'
    return inputlist

def differences(input_list):
    """This function takes the difference between two consecutive items in a list. Also inputs should be float type lists"""
    output_list = []
    for x in range(1,len(input_list)):
        output_list.append(input_list[x]-input_list[x-1])
    return output_list

def eV_to_J(eV):
    return [x*e for x in eV]

def kvalue_adjuster(k_in_list):
    return [x*(2*np.pi/a) for x in k_in_list]

class EffectiveMass():
    """Calculates the mass and does all the calculations"""
    def __init__(self, name):
        self.name = name
        self.file_name_bands = name + ".bands"
        self.file_name_dat = name + ".bands.gnuplot.dat"
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

    def open_file(self):
        """Opens the bands file and reads info"""
        file = open(self.file_name_dat, "r")
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

        #  opening the bands file to get the kpath
        bands_file = open(self.file_name_bands, "r")
        bands_file_rows = []
        hashlinecount = 0  # looking at commented lines in order to extract fermi energy
        for bands_file_line in bands_file:
            bands_file_rows.append(bands_file_line)

        length_bands_file_rows = len(bands_file_rows)
        for i in range(length_bands_file_rows-1, 0, -1):
            temperary = bands_file_rows[i].split()
            if len(temperary) == 2:
                self.kpath_k.append(float(temperary[0].strip()))
                self.kpath_symbol.append(str(temperary[1].strip("'")))
            elif len(temperary) is not 2:
                break
        self.kpath_symbol = latexify(self.kpath_symbol)

    def calculate_mass(self, band_number, minkval, maxkval):
        """This is where the calculation is done"""
        # Identifying the gamma points in the list
        gamma_points = [self.kpath_k[x] for x in range(0,len(self.kpath_symbol)) if self.kpath_symbol[x]=='$\\Gamma$']        
        index = self.band_data_x[band_number].index(gamma_points[0])

        x_to_fit = []
        y_to_fit = []
        xfit = []
        yfit = []
        SecondDeri = []
        effective_masses = []
        for x in gamma_points:
            index = self.band_data_x[band_number].index(x)
            x_temp = []
            y_temp = []
            for datapoint in range(minkval,maxkval+1):
                try:
                    x_temp.append(self.band_data_x[band_number][index + datapoint])
                    y_temp.append(self.band_data_y[band_number][index + datapoint])
                except:
                    pass
            y_temp = eV_to_J(y_temp)
            x_temp = kvalue_adjuster(x_temp)
            x_to_fit.append(x_temp)
            y_to_fit.append(y_temp)

            # Calculate coefficients for the polynomial
            fit_coeffs = np.polyfit(x_temp, y_temp, 2) #making the fit
            SecondDeri.append(fit_coeffs[0])
            effective_masses.append(hbar**2/fit_coeffs[0])
            print(f"Debug: fit_coeffs: {fit_coeffs}")
            print(f"Debug: effective_masses: {effective_masses}")
            print(f"Debug: [x/e for x in effective_masses]: {[x/m0 for x in effective_masses]}")

            #Making the ploynomial
            f = np.poly1d(fit_coeffs)

            # calculate new x's and y's
            x_new = np.linspace(x_temp[0], x_temp[-1], 100)
            y_new = f(x_new)

            xfit.append(x_new)
            yfit.append(y_new)

        plt.figure(1)
        for x in range(0,len(x_to_fit)):
            plt.plot(x_to_fit[x], y_to_fit[x], linewidth=0.4, label=f"Band around $\Gamma$ point")
            plt.plot(xfit[x], yfit[x], linewidth=0.4, label=f"Fit $m^*$: {effective_masses[x]:.5E}")
        plt.xlabel('kpath values')
        plt.ylabel('Energy in J (E - E$_f$)')
        plt.legend(loc='upper right')
        plt.title("Effective mass calculations")
        plt.savefig(f"Band_{band_number}.pdf")

if __name__ == "__main__":
    # Defining all the constants that is used
    Egap = 0.0
    c = 2.99792458E8
    e = 1.60218E-19
    m0 = 9.109E-31
    me = 0.0236*m0
    mh = 0.21*m0
    me_CdSe = 1.18E-31
    mh_CdSe = 4.09E-31
    h = 6.62607015E-34  # (Units Js) Starting from 2019 May 20 This value has been fixed
    hbar = 1.054571800E-34 	
    heV = 4.135667662E-15  # (Units eVs) This is the value is eVs
    epsilon0 = 8.85E-12 # (Units F/m) 
    epsilonr = 24
    a = 6.432E-10 #lattice constant for Sn

    mass = EffectiveMass("Sn")
    mass.open_file()
    mass.calculate_mass(3,-3,3)