import matplotlib.pyplot as plt
import numpy as np

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

def read_scf_out(file_name):
    """
    This function reads the total energy and fermi level from a *.scf.out file.
    |file_name: (str) the file name without the .scf.in part
    |output: float,float
    """
    with open(f"{file_name}.scf.out", "r") as file:
        data = file.read()
        data  = data.split("\n")

        for line in data:
            if "     the Fermi energy is" in line:
                fermi_energy = line.split()[4]
            if "!    total energy              " in line:
                total_energy = line.split()[4]
    return(total_energy, fermi_energy)

class QEBandsReader():
    def __init__(self, file_name):
        print(f"QEBandsReader: Opening Bands file with name '{file_name}'")
        bands_file = open(file_name, "r")
        check = open("check.in", "w+")
        self.data = []
        super_temp = []  # This is where we temporarily store the all k point band data until it is separated into bands
        temp = []  # This is where we temporarily store the individual k point band data until it is separated into bands
        #  Below for loop gets all the per k point band values and saves all of it in the super_temp list
        for row in bands_file:
            # This takes the k points as a list
            if 'number of k points' in row.strip():
                self.k_data = list(range(1,int(row.split()[4])+1))
                # print(self.k_data)
            # This takes the data as a list (we have to invert this list to get actuall plottable data)
            else:
                try:
                    data = row.split()
                    #This initial loop will get rid of any line that has a str 
                    for x in data:
                        float(x)
                    for x in data:
                        temp.append(float(x))
                        # print(temp)
                    # check.write(str(data))
                except:
                    super_temp.append(temp)
                    temp = []
        super_temp = [x for x in super_temp if x != []]
        check.write(str(super_temp))

        check.close()
        # print(len(super_temp[0]))
        # print(super_temp)

        for band_point in range(0,len(super_temp[0])):
            band = []
            for kpoint in range(0,len(super_temp)):
                band.append(super_temp [kpoint][band_point])
            self.data.append(band)

    def calculate_mass(self, band_number, minkval, maxkval):
        """This is where the calculation is done"""
        x_to_fit = []
        y_to_fit = []
        xfit = []
        yfit = []
        SecondDeri = []
        effective_masses = []

        # index_low = self.k_data.index(minkval)
        # index_high = self.k_data.index(maxkval)

        for datapoint in range(minkval,maxkval+1):
            x_to_fit.append(self.k_data[datapoint])
            y_to_fit.append(self.data[band_number][datapoint])

        # Calculate coefficients for the polynomial
        fit_coeffs = np.polyfit(x_to_fit, y_to_fit, 2) #making the fit
        SecondDeri.append(fit_coeffs[0])
        effective_masses.append(hbar**2/fit_coeffs[0])
        print(f"Debug: fit_coeffs: {fit_coeffs}")
        print(f"Debug: effective_masses: {effective_masses}")
        print(f"Debug: [x/e for x in effective_masses]: {[x/m0 for x in effective_masses]}")

        #Making the ploynomial
        f = np.poly1d(fit_coeffs)

        # calculate new x's and y's
        x_new = np.linspace(x_to_fit[0], x_to_fit[-1], 100)
        y_new = f(x_new)

        xfit.append(x_new)
        yfit.append(y_new)

        plt.figure()
        plt.plot(x_to_fit, y_to_fit, linewidth=0.4, label=f"Band around $\\Gamma$ point")
        plt.plot(xfit, yfit, linewidth=0.4,  label=f"Fit $m^*$: {effective_masses[0]:.5E}")
        plt.xlabel('kpath values')
        plt.ylabel('Energy in J (E - E$_f$)')
        plt.legend(loc='upper right')
        plt.title("Effective mass calculations")
        plt.savefig(f"Band_{band_number}.pdf")

if __name__ == "__main__":
    myplot = QEBandsReader("Sn.bands.out")
    # myplot.plot()
    myplot.calculate_mass(14,79,81)