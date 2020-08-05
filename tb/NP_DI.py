"""
This file is meant to be for NPs calulations.
"""

from ebk.tb.TB import TB
import numpy as np
from datetime import datetime
from numpy import linalg as LA
import warnings
warnings.filterwarnings('ignore')

def g(k_point, neighbours):
    g = []
    # Calculating g_0
    val = np.exp(np.dot(k_point,neighbours[0])*1j) + np.exp(np.dot(k_point,neighbours[1])*1j) + np.exp(np.dot(k_point,neighbours[2])*1j) + np.exp(np.dot(k_point,neighbours[3])*1j)
    g.append(val)
    # Calculating g_1
    val = np.exp(np.dot(k_point,neighbours[0])*1j) + np.exp(np.dot(k_point,neighbours[1])*1j) - np.exp(np.dot(k_point,neighbours[2])*1j) - np.exp(np.dot(k_point,neighbours[3])*1j)
    g.append(val)
    # Calculating g_2
    val = np.exp(np.dot(k_point,neighbours[0])*1j) - np.exp(np.dot(k_point,neighbours[1])*1j) + np.exp(np.dot(k_point,neighbours[2])*1j) - np.exp(np.dot(k_point,neighbours[3])*1j)
    g.append(val)
    # Calculating g_3
    val = np.exp(np.dot(k_point,neighbours[0])*1j) - np.exp(np.dot(k_point,neighbours[1])*1j) - np.exp(np.dot(k_point,neighbours[2])*1j) + np.exp(np.dot(k_point,neighbours[3])*1j)
    g.append(val)
    return g

def distance(v1, v2):
    """Gives the distance vectors for v1-v2"""
    v = []
    for i in range(1,4):
        v.append(v1[i]-v2[i])
    return v, np.sqrt(v[0]**2+v[1]**2+v[2]**2)

class DI(TB):
    """
    This class is for nanoparticles which are ZB like and has an extra atom (like H) for passivation.
    Here naturally we do a single k point calculation.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def plot(self):
        import matplotlib.pyplot as plt
        x = [0,1]
        x_0 = [0,1.1]
        E_HOMO = self.sorted_eigen_vals[self.HOMO_i]
        for i, val in enumerate(self.sorted_eigen_vals):
            y = [val, val]
            if i == self.HOMO_i:
                print(f"homo level {val} which is the {i}th index")
                plt.plot(x_0, y, "r", label = "HOMO level")
            else:
                plt.plot(x, y, "b")
        plt.title("Energies")
        plt.ylabel = f"Energy (eV)"
        plt.xticks(x, " ")
        plt.legend()
        plt.savefig(f"{self.run_name}.pdf")
        plt.ylim(E_HOMO - 3, E_HOMO + 3)
        # plt.show()
        plt.savefig(f"{self.run_name}_zoomed.pdf")

    def load_parameters(self):
        """This method will eventually be able to load from file. Since we dont need it now we can just do a quick assignment of the parameters"""
        # We save all the parameters here into a dictionary
        self.parameter = {}

        if self.so == True:
            self.so_a = 0.25  # This is in eV
            self.so_c = 0.25  # This is in eV
        elif self.so == False:
            self.so_a = 0.0  # This is in eV
            self.so_c = 0.0  # This is in eV

        # This is for HgTe
        # self.parameter.update({"so_a": 2*self.so_a})
        self.parameter.update({"so_a": 2*self.so_a})
        self.parameter.update({"so_c": self.so_c})

        # Defining tight bindging paramters here _ca only for diagonal elements (remember we are using the ac block matrix as the base and then getting ca just by complex conjugate transpose)
        # The anion energies
        self.parameter.update({"Es_a": -7.2029})
        self.parameter.update({"Ep_a": 0.0032})
        self.parameter.update({"Edxy_a": 0})
        self.parameter.update({"Edx2my2_a":0})
        self.parameter.update({"ES_a": 14.1081})

        # The cation energies
        # self.parameter.update({"Es_c": self.parameter["Es_a"]})
        # self.parameter.update({"Ep_c": self.parameter["Ep_a"]})
        # self.parameter.update({"Edxy_c": self.parameter["Edxy_a"]})
        # self.parameter.update({"Edx2my2_c": self.parameter["Edx2my2_a"]})
        # self.parameter.update({"ES_c": self.parameter["ES_a"]})

        # Now we load the other parameters
        # The diagonal parameters
        self.parameter.update({"Vsssigma": -1.4802})
        self.parameter.update({"Vddsigma": 0})  # meka thamai missing thiyenne 
        self.parameter.update({"Vppsigma": 3.3555})
        self.parameter.update({"Vpppi": -1.6848})
        self.parameter.update({"Vddpi": 0})
        self.parameter.update({"Vdddelta": 0})
        self.parameter.update({"VSSsigma": 0.4795})

        # The anion to cation off-diagonal parameters
        self.parameter.update({"Vspsigma_ac": 1.9034})
        self.parameter.update({"Vsdsigma_ac": 0})
        self.parameter.update({"Vpdsigma_ac": 0})
        self.parameter.update({"Vpdpi_ac": 0})
        self.parameter.update({"VsS_ac": 0})
        self.parameter.update({"VSpsigma_ac": 2.7564})
        self.parameter.update({"VSdsigma_ac":0})

        # # The cation to anion off-diagonal parameters
        # self.parameter.update({"Vspsigma_ca": 1.79243})
        # self.parameter.update({"Vsdsigma_ca": -2.06445})
        # self.parameter.update({"Vpdsigma_ca": -1.27510})
        # self.parameter.update({"Vpdpi_ca": 1.35105})
        # self.parameter.update({"VsS_ca": 0.62601})
        # self.parameter.update({"VSpsigma_ca": 0.62509})
        # self.parameter.update({"VSdsigma_ca": 0.73680})

        # The Hydrogen energies with itself
        self.parameter.update({"E_H": 0})

        # The Hydrogen diagonal energies
        self.parameter.update({"Vsssigma_H": 0})

        # The Hydrogen off diagonal energies
        self.parameter.update({"Vspsigma_H": 0})



        # if self.so == True:
        #     self.so_a = 0.25  # This is in eV
        #     self.so_c = 0.25  # This is in eV
        # elif self.so == False:
        #     self.so_a = 0.0  # This is in eV
        #     self.so_c = 0.0  # This is in eV

        # # This is for HgTe
        # # self.parameter.update({"so_a": 2*self.so_a})
        # self.parameter.update({"so_a": self.so_a})
        # self.parameter.update({"so_c": self.so_c})

        # # Defining tight bindging paramters here _ca only for diagonal elements (remember we are using the ac block matrix as the base and then getting ca just by complex conjugate transpose)
        # # The anion energies
        # self.parameter.update({"Es_a": -5.32567})
        # self.parameter.update({"Ep_a": 2.52194})
        # self.parameter.update({"Edxy_a": 12.06317})
        # self.parameter.update({"Edx2my2_a": 12.06317})
        # self.parameter.update({"ES_a": 8.80303})

        # # The cation energies
        # self.parameter.update({"Es_c": -5.32567})
        # self.parameter.update({"Ep_c": 2.52194})
        # self.parameter.update({"Edxy_c": 12.06317})
        # self.parameter.update({"Edx2my2_c": 12.06317})
        # self.parameter.update({"ES_c": 8.80303})

        # # Now we load the other parameters
        # # The diagonal parameters
        # self.parameter.update({"Vsssigma": -1.26716})
        # self.parameter.update({"Vddsigma": -2.33598})  # meka thamai missing thiyenne 
        # self.parameter.update({"Vppsigma": 2.75609})
        # self.parameter.update({"Vpppi": -1.11032})
        # self.parameter.update({"Vddpi": 2.53095})
        # self.parameter.update({"Vdddelta": -1.85318})
        # self.parameter.update({"VSSsigma": -0.93471})

        # # The anion to cation off-diagonal parameters
        # self.parameter.update({"Vspsigma_ac": 1.79243})
        # self.parameter.update({"Vsdsigma_ac": -2.06445})
        # self.parameter.update({"Vpdsigma_ac": -1.27510})
        # self.parameter.update({"Vpdpi_ac": 1.35105})
        # self.parameter.update({"VsS_ac": 0.62601})
        # self.parameter.update({"VSpsigma_ac": 0.62509})
        # self.parameter.update({"VSdsigma_ac":0.73680})

        # # The cation to anion off-diagonal parameters
        # self.parameter.update({"Vspsigma_ca": 1.79243})
        # self.parameter.update({"Vsdsigma_ca": -2.06445})
        # self.parameter.update({"Vpdsigma_ca": -1.27510})
        # self.parameter.update({"Vpdpi_ca": 1.35105})
        # self.parameter.update({"VsS_ca": 0.62601})
        # self.parameter.update({"VSpsigma_ca": 0.62509})
        # self.parameter.update({"VSdsigma_ca": 0.73680})

        # # The Hydrogen energies with itself
        # self.parameter.update({"E_H": 0.56000})

        # # The Hydrogen diagonal energies
        # self.parameter.update({"Vsssigma_H": -5.15100})

        # # The Hydrogen off diagonal energies
        # self.parameter.update({"Vspsigma_H": 5.27000})

        # Save new list with all paramters converted to J
        # self.parameterJ = {key: value * 1.60218e-19 for (key,value) in self.parameter.items()}
        self.parameterJ = {key: value for (key,value) in self.parameter.items()}

    def Vsp(self, direction):
        if direction == "ac":
            return self.parameterJ["Vspsigma_ac"]/np.sqrt(3)
        elif direction == "ca":
            return self.parameterJ["Vspsigma_ca"]/np.sqrt(3)

    def Vhp(self, direction):
        if direction == "ac":
            return self.parameterJ["Vspsigma_H"]/np.sqrt(3)
        elif direction == "ca":
            # This part is not implemented no harm in leaving it here since we dont really use it.
            return self.parameterJ["Vspsigma_ca"]/np.sqrt(3)

    def Vsd(self, direction):
        if direction == "ac":
            Vpdsigma = self.parameterJ["Vsdsigma_ac"]
        elif direction == "ca":
            Vpdsigma = self.parameterJ["Vsdsigma_ca"]
        return Vpdsigma/np.sqrt(3)

    def Vpd(self, direction):
        if direction == "ac":
            Vpdsigma = self.parameterJ["Vpdsigma_ac"]
            Vpdpi = self.parameterJ["Vpdpi_ac"]
        elif direction == "ca":
            Vpdsigma = self.parameterJ["Vpdsigma_ca"]
            Vpdpi = self.parameterJ["Vpdpi_ca"]
        return (np.sqrt(3)*Vpdsigma - 2*Vpdpi)/(3*np.sqrt(3))

    def Vxx(self):
        Vppsigma = self.parameterJ["Vppsigma"]
        Vpppi = self.parameterJ["Vpppi"]
        return (1/3)*Vppsigma + (2/3)*Vpppi

    def Vxy(self):
        Vppsigma = self.parameterJ["Vppsigma"]
        Vpppi = self.parameterJ["Vpppi"]
        return (Vppsigma - Vpppi)/3

    def Upd(self, direction):
        if direction == "ac":
            Vpdsigma = self.parameterJ["Vpdsigma_ac"]
            Vpdpi = self.parameterJ["Vpdpi_ac"]
        elif direction == "ca":
            Vpdsigma = self.parameterJ["Vpdsigma_ca"]
            Vpdpi = self.parameterJ["Vpdpi_ca"]
        return (1/(3*np.sqrt(3)))*(np.sqrt(3)*Vpdsigma + Vpdpi)

    def Wpd(self, direction):
        if direction == "ac":
            Vpdpi = self.parameterJ["Vpdpi_ac"]
        elif direction == "ca":
            Vpdpi = self.parameterJ["Vpdpi_ca"]
        return 2*(Vpdpi)/3

    def Vdd(self):
        Vddsigma = self.parameterJ["Vddsigma"]
        Vddpi = self.parameterJ["Vddpi"]
        Vdddelta = self.parameterJ["Vdddelta"]
        return (1/3)*(Vddsigma + (2/3)*Vddpi + (4/3)*Vdddelta)

    def Vdd_2(self):
        Vddpi = self.parameterJ["Vddpi"]
        Vdddelta = self.parameterJ["Vdddelta"]
        return (1/3)*(2*Vddpi + Vdddelta)

    def Udd(self):
        Vddsigma = self.parameterJ["Vddsigma"]
        Vddpi = self.parameterJ["Vddpi"]
        Vdddelta = self.parameterJ["Vdddelta"]
        return (1/3)*(Vddsigma - Vddpi - (2/3)*Vdddelta)

    def Wdd(self):
        Vddpi = self.parameterJ["Vddpi"]
        Vdddelta = self.parameterJ["Vdddelta"]
        return (2/(3*np.sqrt(3)))*(Vdddelta - Vddpi)

    def VSd(self, direction):
        if direction == "ac":
            VSd = self.parameterJ["VSdsigma_ac"]
        elif direction == "ca":
            VSd = self.parameterJ["VSdsigma_ca"]
        return VSd/np.sqrt(3)

    def VSp(self, direction):
        if direction == "ac":
            VSp = self.parameterJ["VSpsigma_ac"]
        elif direction == "ca":
            VSp = self.parameterJ["VSpsigma_ca"]
        return VSp/np.sqrt(3)

    def VsS(self, direction):
        if direction == "ac":
            VsS = self.parameterJ["VsS_ac"]
        elif direction == "ca":
            VsS = self.parameterJ["VsS_ca"]
        return VsS

    def function_chooser(self, row, col, k_point = [0, 0, 0]):
        if row == col:
            # print(row, col)
            # Here are the diagonal elements of the big matrix
            if row <= 9:
            # These are the anion to anion elements
            if row == 0: return self.parameterJ["Es_a"]
            elif row == 1 or row == 2 or row ==3 : return self.parameterJ["Ep_a"]
            elif row == 4 or row == 5 or row ==6: return self.parameterJ["Edxy_a"]
            elif row == 7 or row == 8: return self.parameterJ["Edx2my2_a"]
            elif row == 9: return self.parameterJ["ES_a"]
            else:
                print("function_chooser: Fatal Error - Element not found")

            elif row > 9:
                row = row - 10  # So we can still refer to the rows and columns as from 0-9 to retain our sanity
                # These are the cation to cation elements
                if row == 0: return self.parameterJ["Es_c"]
                elif row == 1 or row == 2 or row ==3 : return self.parameterJ["Ep_c"]
                elif row == 4 or row == 5 or row ==6: return self.parameterJ["Edxy_c"]
                elif row == 7 or row == 8: return self.parameterJ["Edx2my2_c"]
                elif row == 9: return self.parameterJ["ES_c"]
        elif row <= 9 and col > 9:
            # These are the anion to cation elements
            col = col - 10  # So we can still refer to the rows and columns as from 0-9 to retain our sanity
            # next let's select the diagonal and upper triangle
            neighbours_ac = self.cations_as_neighbours  # this has to be edited to make it general and have to consider the two atoms.

            if row == 0:
                if col == 0: return self.parameterJ["Vsssigma"]*g(k_point,neighbours_ac)[0]
                elif col == 1: return self.Vsp("ac")*g(k_point,neighbours_ac)[1]
                elif col == 2: return self.Vsp("ac")*g(k_point,neighbours_ac)[2]
                elif col == 3: return self.Vsp("ac")*g(k_point,neighbours_ac)[3]
                elif col == 4: return self.Vsd("ac")*g(k_point,neighbours_ac)[3]
                elif col == 5: return self.Vsd("ac")*g(k_point,neighbours_ac)[1]
                elif col == 6: return self.Vsd("ac")*g(k_point,neighbours_ac)[2]
                elif col == 9: return self.VsS("ac")*g(k_point,neighbours_ac)[0]
                else:
                    return 0. + 0.j

            elif row == 1:
                if col == 0: return -self.Vsp("ca")*g(k_point,neighbours_ac)[1]
                elif col == 1: return self.Vxx()*g(k_point,neighbours_ac)[0]
                elif col == 2: return self.Vxy()*g(k_point,neighbours_ac)[3]
                elif col == 3: return self.Vxy()*g(k_point,neighbours_ac)[2]
                elif col == 4: return self.Upd("ac")*g(k_point,neighbours_ac)[2]
                elif col == 5: return self.Vpd("ac")*g(k_point,neighbours_ac)[0]
                elif col == 6: return self.Upd("ac")*g(k_point,neighbours_ac)[3]
                elif col == 7: return (np.sqrt(3)/2)*self.Wpd("ac")*g(k_point,neighbours_ac)[1]
                elif col == 8: return -(1/2)*self.Wpd("ac")*g(k_point,neighbours_ac)[1]
                elif col == 9: return -self.VSp("ca")*g(k_point,neighbours_ac)[1]
                else:
                    return 0. + 0.j

            elif row == 2:
                if col == 0: return -self.Vsp("ca")*g(k_point,neighbours_ac)[2]
                elif col == 1: return self.Vxy()*g(k_point,neighbours_ac)[3]
                elif col == 2: return self.Vxx()*g(k_point,neighbours_ac)[0]
                elif col == 3: return self.Vxy()*g(k_point,neighbours_ac)[1]
                elif col == 4: return self.Upd("ac")*g(k_point,neighbours_ac)[1]
                elif col == 5: return self.Upd("ac")*g(k_point,neighbours_ac)[3]
                elif col == 6: return self.Vpd("ac")*g(k_point,neighbours_ac)[0]
                elif col == 7: return -(np.sqrt(3)/2)*self.Wpd("ac")*g(k_point,neighbours_ac)[2]
                elif col == 8: return -(1/2)*self.Wpd("ac")*g(k_point,neighbours_ac)[2]
                elif col == 9: return -self.VSp("ca")*g(k_point,neighbours_ac)[2]
                else:
                    return 0. + 0.j

            elif row == 3:
                if col == 0: return -self.Vsp("ca")*g(k_point,neighbours_ac)[3]
                elif col == 1: return self.Vxy()*g(k_point,neighbours_ac)[2]
                elif col == 2: return self.Vxy()*g(k_point,neighbours_ac)[1]
                elif col == 3: return self.Vxx()*g(k_point,neighbours_ac)[0]
                elif col == 4: return self.Vpd("ac")*g(k_point,neighbours_ac)[0]
                elif col == 5: return self.Upd("ac")*g(k_point,neighbours_ac)[2]
                elif col == 6: return self.Upd("ac")*g(k_point,neighbours_ac)[1]
                elif col == 8: return self.Wpd("ac")*g(k_point,neighbours_ac)[3]
                elif col == 9: return -self.VSp("ca")*g(k_point,neighbours_ac)[3]
                else:
                    return 0. + 0.j

            if row == 4:
                if col == 0: return self.Vsd("ca")*g(k_point,neighbours_ac)[3]
                elif col == 1: return -self.Upd("ca")*g(k_point,neighbours_ac)[2]
                elif col == 2: return -self.Upd("ca")*g(k_point,neighbours_ac)[1]
                elif col == 3: return -self.Vpd("ca")*g(k_point,neighbours_ac)[0]
                elif col == 4: return self.Vdd()*g(k_point,neighbours_ac)[0]
                elif col == 5: return self.Udd()*g(k_point,neighbours_ac)[2]
                elif col == 6: return self.Udd()*g(k_point,neighbours_ac)[1]
                elif col == 8: return self.Wdd()*g(k_point,neighbours_ac)[3]
                elif col == 9: return self.VSd("ca")*g(k_point,neighbours_ac)[3]
                else:
                    return 0. + 0.j

            if row == 5:
                if col == 0: return self.Vsd("ca")*g(k_point,neighbours_ac)[1]
                elif col == 1: return -self.Vpd("ca")*g(k_point,neighbours_ac)[0]
                elif col == 2: return -self.Upd("ca")*g(k_point,neighbours_ac)[3]
                elif col == 3: return -self.Upd("ca")*g(k_point,neighbours_ac)[2]
                elif col == 4: return self.Udd()*g(k_point,neighbours_ac)[2]
                elif col == 5: return self.Vdd()*g(k_point,neighbours_ac)[0]
                elif col == 6: return self.Udd()*g(k_point,neighbours_ac)[3]
                elif col == 7: return (np.sqrt(3)/2)*self.Wdd()*g(k_point,neighbours_ac)[1]
                elif col == 8: return -(1/2)*self.Wdd()*g(k_point,neighbours_ac)[1]
                elif col == 9: return self.VSd("ca")*g(k_point,neighbours_ac)[1]
                else:
                    return 0. + 0.j

            if row == 6:
                if col == 0: return self.Vsd("ca")*g(k_point,neighbours_ac)[2]
                elif col == 1: return -self.Upd("ca")*g(k_point,neighbours_ac)[3]
                elif col == 2: return -self.Vpd("ca")*g(k_point,neighbours_ac)[0]
                elif col == 3: return -self.Upd("ca")*g(k_point,neighbours_ac)[1]
                elif col == 4: return self.Udd()*g(k_point,neighbours_ac)[1]
                elif col == 5: return self.Udd()*g(k_point,neighbours_ac)[3]
                elif col == 6: return self.Vdd()*g(k_point,neighbours_ac)[0]
                elif col == 7: return -(np.sqrt(3)/2)*self.Wdd()*g(k_point,neighbours_ac)[2]
                elif col == 8: return -(1/2)*self.Wdd()*g(k_point,neighbours_ac)[2]
                elif col == 9: return self.VSd("ca")*g(k_point,neighbours_ac)[2]
                else:
                    return 0. + 0.j

            if row == 7:
                if col == 1: return -(np.sqrt(3)/2)*self.Wpd("ca")*g(k_point,neighbours_ac)[1]
                elif col == 2: return (np.sqrt(3)/2)*self.Wpd("ca")*g(k_point,neighbours_ac)[2]
                elif col == 5: return (np.sqrt(3)/2)*self.Wdd()*g(k_point,neighbours_ac)[1]
                elif col == 6: return -(np.sqrt(3)/2)*self.Wdd()*g(k_point,neighbours_ac)[2]
                elif col == 7: return self.Vdd_2()*g(k_point,neighbours_ac)[0]
                else:
                    return 0. + 0.j

            if row == 8:
                if col == 1: return (1/2)*self.Wpd("ca")*g(k_point,neighbours_ac)[1]
                elif col == 2: return (1/2)*self.Wpd("ca")*g(k_point,neighbours_ac)[2]
                elif col == 3: return -self.Wpd("ca")*g(k_point,neighbours_ac)[3]
                elif col == 4: return self.Wdd()*g(k_point,neighbours_ac)[3]
                elif col == 5: return -(1/2)*self.Wdd()*g(k_point,neighbours_ac)[1]
                elif col == 6: return -(1/2)*self.Wdd()*g(k_point,neighbours_ac)[2]
                elif col == 8: return self.Vdd_2()*g(k_point,neighbours_ac)[0]
                else:
                    return 0. + 0.j

            if row == 9:
                if col == 0: return self.VsS("ca")*g(k_point,neighbours_ac)[0]
                elif col == 1: return self.VSp("ac")*g(k_point,neighbours_ac)[1]
                elif col == 2: return self.VSp("ac")*g(k_point,neighbours_ac)[2]
                elif col == 3: return self.VSp("ac")*g(k_point,neighbours_ac)[3]
                elif col == 4: return self.VSd("ac")*g(k_point,neighbours_ac)[3]
                elif col == 5: return self.VSd("ac")*g(k_point,neighbours_ac)[1]
                elif col == 6: return self.VSd("ac")*g(k_point,neighbours_ac)[2]
                elif col == 9: return self.parameterJ["VSSsigma"]*g(k_point,neighbours_ac)[0]
                else:
                    return 0. + 0.j
            else:
                return 0. + 0.j

        elif row > 9 and col <= 9:
            # this part is not implemented in our calculationas and is done elsewhere.
            # These are the cations to anion elements
            # They will be the conjugate transpose of the anion to cation elements
            # Therefore we have to make sure that this will by done after the row<=9 and col > 9 set is done with
            return np.conjugate(self.H[col][row])
        else:
            return 0. + 0.j

    def Haa_f(self):
        """This is defined for the anion"""
        self.Haa = np.zeros((10, 10), dtype=complex)
        for i in range(10):
            self.Haa[i][i] = self.function_chooser(i,i)
        self.log("----- Haa block:")
        self.log(f"{self.Haa}")

    def Hcc_f(self):
        """This is defined for the cation"""
        self.Hcc = np.zeros((10, 10), dtype=complex)
        for i in range(10,20):
            self.Hcc[i-10][i-10] = self.function_chooser(i,i)
        self.log("----- Hcc block:")
        self.log(f"{self.Hcc}")

    def Hac_f(self):
        """This is defined for the anion to cation elements"""
        self.Hac = np.zeros((10, 10), dtype=complex)
        for i in range(10):
            for j in range(10):
                self.Hac[i][j] = self.function_chooser(i,j+10)
        self.log("----- Hac block:")
        self.log(f"{self.Hac}")

    def Hca_f(self):
        """This is defined for the anion to cation elements
        To implement this matrix Hac should already be caclulated"""
        self.Hca = np.zeros((10, 10), dtype=complex)
        for i in range(10):
            for j in range(10):
                self.Hca[i][j] = np.conjugate(self.Hac[j][i])
        self.log("----- Hca block:")
        self.log(f"{self.Hca}")

    def Hhh_f(self):
        self.Hhh = np.zeros((1, 1), dtype=complex)
        for i in range(1):
            self.Hhh[i][i] = self.parameterJ["E_H"]
        self.log("----- Hhh block:")
        self.log(f"{self.Hhh}")
        
    def Hha_f(self):
        self.Hha = np.zeros((1, 10), dtype=complex)
        for i in range(1):
            for col in range(10):
                if col == 0: self.Hha[i][col] = self.parameterJ["Vsssigma_H"]*g([0,0,0],self.cations_as_neighbours)[0]
                elif col == 1: self.Hha[i][col] = self.Vhp("ac")*g([0, 0, 0],self.cations_as_neighbours)[1]
                elif col == 2: self.Hha[i][col] = self.Vhp("ac")*g([0, 0, 0],self.cations_as_neighbours)[2]
                elif col == 3: self.Hha[i][col] = self.Vhp("ac")*g([0, 0, 0],self.cations_as_neighbours)[3]
                else: self.Hha[i][col] = 0. + 0.j
        self.log("----- Hha block:")
        self.log(f"{self.Hha}")

    def Hah_f(self):
        self.Hah = np.zeros((10, 1), dtype=complex)
        for i in range(10):
            for j in range(1):
                self.Hah[i][j] = np.conjugate(self.Hha[j][i])
        self.log("----- Hah block:")
        self.log(f"{self.Hah}")

    def Hch_f(self):
        self.Hch = self.Hah.copy()
        self.log("----- Hch block:")
        self.log(f"{self.Hch}")

    def Hhc_f(self):
        self.Hhc = self.Hha.copy()
        self.log("----- Hhc block:")
        self.log(f"{self.Hhc}\n")

    def get_r_c_d(self, s_r, s_c , n_r, n_c, o_r, o_c):
        """This function outputs the columns and rows and lso the relevant SO term for that element
        you still have to feed in the correct o_r and o_c depending whether it is a H or not
        """
        D = 0
        if n_r == n_c:
            # Spin orbit couling is for the same atom.
            # print(f"{n_r} and {n_c}, sr:{s_r} sc:{s_c}, oc:{o_c} or:{o_r}")
            # print(n_r == n_c)
            if s_r == 0 and s_c == 0:
                # print(f"{n_r} and {n_c}, sr:{s_r} sc:{s_c}, oc:{o_c} or:{o_r}")
                # print(f"inside loop")
                # Doing the up up submatrix
                if self.atoms[n_c][0] == "Sn":
                    # Doing the anions
                    if o_r == 2 and o_c == 1: D = self.parameterJ["so_a"]*1j # py, px
                    if o_r == 1 and o_c == 2: D = -self.parameterJ["so_a"]*1j # px, py
                # self.SuperH[self.i_up+self.i_anion+self.i_py][self.i_up+self.i_anion+self.i_px] += self.parameterJ["so_a"]*1j
                # self.SuperH[self.i_up+self.i_anion+self.i_px][self.i_up+self.i_anion+self.i_py] += -self.parameterJ["so_a"]*1j
                if self.atoms[n_c][0] == "Sn":
                    # Doing the cations
                    if o_r == 2 and o_c == 1: D = self.parameterJ["so_c"]*1j  # py, px
                    if o_r == 1 and o_c == 2: D = -self.parameterJ["so_c"]*1j # px, py
                # self.SuperH[self.i_up+self.i_cation+self.i_py][self.i_up+self.i_cation+self.i_px] += self.parameterJ["so_c"]*1j
                # self.SuperH[self.i_up+self.i_cation+self.i_px][self.i_up+self.i_cation+self.i_py] += -self.parameterJ["so_c"]*1j
            if s_r == 1 and s_c == 1:
                # Doing the down down submatrix
                if self.atoms[n_c][0] == "Sn":
                    # Doing the anions
                    if o_r == 2 and o_c == 1: D = -self.parameterJ["so_a"]*1j  # py, px
                    if o_r == 1 and o_c == 2: D = self.parameterJ["so_a"]*1j  # px, py
                # # Lets do the down-down sub-matrix next
                # self.SuperH[self.i_down+self.i_anion+self.i_py][self.i_down+self.i_anion+self.i_px] += -self.parameterJ["so_a"]*1j
                # self.SuperH[self.i_down+self.i_anion+self.i_px][self.i_down+self.i_anion+self.i_py] += self.parameterJ["so_a"]*1j
                if self.atoms[n_c][0] == "Sn":
                    # Doing the cations
                    if o_r == 2 and o_c == 1: D = -self.parameterJ["so_c"]*1j  # py, px
                    if o_r == 1 and o_c == 2: D = self.parameterJ["so_c"]*1j # px, py
                # self.SuperH[self.i_down+self.i_cation+self.i_py][self.i_down+self.i_cation+self.i_px] += -self.parameterJ["so_c"]*1j
                # self.SuperH[self.i_down+self.i_cation+self.i_px][self.i_down+self.i_cation+self.i_py] += self.parameterJ["so_c"]*1j
            if s_r == 0 and s_c == 1:
                # Doing the up down submatrix
                if self.atoms[n_c][0] == "Sn":
                    # Doing the anions
                    if o_r == 3 and o_c == 1: D = -self.parameterJ["so_a"]
                    if o_r == 1 and o_c == 3: D = self.parameterJ["so_a"]
                    if o_r == 3 and o_c == 2: D = self.parameterJ["so_a"]*1j
                    if o_r == 2 and o_c == 3: D = -self.parameterJ["so_a"]*1j
                if self.atoms[n_c][0] == "Sn":
                    # Doing the cations
                    if o_r == 3 and o_c == 1: D = -self.parameterJ["so_c"]
                    if o_r == 1 and o_c == 3: D = self.parameterJ["so_c"]
                    if o_r == 3 and o_c == 2: D = self.parameterJ["so_c"]*1j
                    if o_r == 2 and o_c == 3: D = -self.parameterJ["so_c"]*1j
            if s_r == 1 and s_c == 0:
                # Doing the down up submatrix
                if self.atoms[n_c][0] == "Sn":
                    # Doing the anions
                    if o_r == 1 and o_c == 3: D = -self.parameterJ["so_a"]
                    if o_r == 3 and o_c == 1: D = self.parameterJ["so_a"]
                    if o_r == 2 and o_c == 3: D = -self.parameterJ["so_a"]*1j
                    if o_r == 3 and o_c == 2: D = self.parameterJ["so_a"]*1j
                if self.atoms[n_c][0] == "Sn":
                    # Doing the cations
                    if o_r == 1 and o_c == 3: D = -self.parameterJ["so_c"]
                    if o_r == 3 and o_c == 1: D = self.parameterJ["so_c"]
                    if o_r == 2 and o_c == 3: D = -self.parameterJ["so_c"]*1j
                    if o_r == 3 and o_c == 2: D = self.parameterJ["so_c"]*1j

        row = 0
        for i in range(n_r):
            # for j in range(s_r+1):
            if self.atoms[i][0] == "Sn" or self.atoms[i][0] == "Sn":
                row+=10
            if self.atoms[i][0] == "H":
                row+=1

        col = 0
        for i in range(n_c):
            # for j in range(s_c+1):
            if self.atoms[i][0] == "Sn" or self.atoms[i][0] == "Sn":
                col+=10
            if self.atoms[i][0] == "H":
                col+=1

        if s_r == 1: row+= int(self.size/2)
        if s_c == 1: col+= int(self.size/2)


        test_row = row + o_r
        test_col = col + o_c
        if (test_row == 0 and test_col == 9) or (test_row == 9 and test_col == 0):
            print(s_r, s_c, n_r, n_c, o_r, o_c)

        return row + o_r, col + o_c, D

    def calculate(self):
        if self.check():
            self.log("**** Calculation has not finished.....")
            self.log(f"**** Possible errors: {self.check}")
            self.log("**** Refer to out file to see error")

        # Here we calculate teh shape of the matrix
        # Te is taken to be the anion in the paper and Hg is taken to be the cation
        # Counting the number of anions and cations and mutiplying by 10 since they all have spds* (10 orbtal) orbitals
        self.anion_count = 0
        self.cation_count = 0
        self.H_count = 0

        for atom in self.atoms:
            if atom[0] == "Sn": self.anion_count+=1
            elif atom[0] == "H": self.H_count+=1
            
        
        self.log("--------- Atom counts:")
        self.log(f"          Dot  (Sn)\t\t:{self.anion_count}\t-> electrons: 6*{self.anion_count} = {6*self.anion_count}") 
        self.log(f"          Hydrogens\t\t:{self.H_count}")
        self.log(f"Size of matrix\t: 2*[(anions:{self.anion_count})*10 + (H atoms:{self.H_count})*1]\n")

        size_wo_so = (self.anion_count*10 + self.H_count)
        size_so = 2*(self.anion_count*10 + self.H_count)

        if self.so == True:
            self.size = size_so
            self.log("This is a spin_orbit enabled calculation.")
        else:
            self.size = size_wo_so

        self.SuperH = np.zeros((self.size, self.size), dtype=complex)
        self.log(f"Shape of initialized matrix : {self.SuperH.shape}\n")
        self.HOMO_i = self.anion_count*6+self.cation_count*2+self.H_count*1 - 1  # -1 for the index
        self.log(f"HOMO level at {self.HOMO_i + 1}th Electronic state\n")

        # print(self.SuperH)
        # calculating priliminary matrices that will be used.
        self.load_parameters()
        self.Haa_f()
        self.Hcc_f()
        self.Hac_f()
        self.Hca_f()
        self.Hhh_f()
        self.Hha_f()
        self.Hah_f()
        self.Hch_f()
        self.Hhc_f()
        print("Done with premliminary matrices")

        self.log("--------- Printing Loaded parameters (J):")
        for i, (k,v) in enumerate(self.parameterJ.items()):
            self.log(f"{k:<11}: {v:>10}")
        self.log("--------- End of Loaded parameters")

        self.log(f"\nNearest neighbour threshold distance: {self.nn_threshold}")
        # self.log("--------- Nearest neighbours:")

        # Building the matrix
        print(f"Building the super matrix")
        for s_r in self.S:  # Spin_rows
            for s_c in self.S:  # Spin_columns
                for n_r in self.N:  # iterating over the atoms for rows
                    for n_c in self.N:  # iterating over atoms for columns
                # first checking if they are nearest neighbours
                        if s_c >= s_r and n_c >= n_r:
                            v, d = distance(self.atoms[n_r],self.atoms[n_c])
                            if d <= self.nn_threshold:
                                if n_r == n_c:
                                    # Case where we have to use Haa, Hcc or Hhh matrices
                                    if self.atoms[n_r][0] == "Sn" and self.atoms[n_c][0] == "Sn":  # Making sure it is not a H
                                        # print(f"n_r is: {n_r} and s_r is {s_r} and s_c is {s_c}")
                                        for o_r in self.O:
                                            for o_c in self.O:
                                                i, j, D = self.get_r_c_d(s_r, s_c , n_r, n_c, o_r, o_c)
                                                # print(i,j)
                                                if self.atoms[n_r][0] == "Sn":
                                                    if s_c == s_r: self.SuperH[i][j] = self.Haa[o_r][o_c] + D
                                                    else: self.SuperH[i][j] = D
                                                    self.SuperH[j][i] = np.conjugate(self.SuperH[i][j])
                                                    # if s_r == 1 and s_c == 1:
                                                        # print(f"{self.Hcc[o_r][o_c]}")
                                                # if self.atoms[n_r][0] == "Sn":
                                                #     if s_c == s_r: self.SuperH[i][j] = self.Haa[o_r][o_c] + D
                                                #     else: self.SuperH[i][j] = D
                                                # self.SuperH[j][i] = np.conjugate(self.SuperH[i][j])
                                    if self.atoms[n_r][0] == "H":
                                        # print(f"n_r is: {n_r} and s_r is {s_r} and s_c is {s_c}")
                                        for o_r in self.H:
                                            for o_c in self.H:
                                                i, j, D = self.get_r_c_d(s_r, s_c , n_r, n_c, o_r, o_c)
                                                if s_c == s_r: self.SuperH[i][j] = self.Hhh[o_r][o_c] + D
                                                else: self.SuperH[i][j] = D
                                                self.SuperH[j][i] = np.conjugate(self.SuperH[i][j])
                                    # self.log(f"nns: {self.atoms[n_r]} and {self.atoms[n_c]}")

                                if self.atoms[n_r][0] == "Sn":
                                    # This is a Hac type block
                                    if self.atoms[n_c][0] == "Sn":  # Hg is the cation
                                        # print(f"n_r is: {n_r} and s_r is {s_r} and s_c is {s_c}")
                                        for o_r in self.O:
                                            for o_c in self.O:
                                                # Here we assume that since we are considering only nearest neighbours we have the other type of atoms in teh coulmns
                                                i, j, D = self.get_r_c_d(s_r, s_c , n_r, n_c, o_r, o_c)
                                                if s_c == s_r: self.SuperH[i][j] = self.Haa[o_r][o_c] + D
                                                else: self.SuperH[i][j] = D
                                                self.SuperH[j][i] = np.conjugate(self.SuperH[i][j])
                                        # self.log(f"nns: {self.atoms[n_r]} and {self.atoms[n_c]}")
                                    if self.atoms[n_c][0] == "H":
                                        # print(f"n_r is: {n_r} and s_r is {s_r} and s_c is {s_c}")
                                        for o_r in self.O:
                                            for o_c in self.H:
                                                i, j, D = self.get_r_c_d(s_r, s_c , n_r, n_c, o_r, o_c)
                                                if s_c == s_r: self.SuperH[i][j] = self.Hah[o_r][o_c] + D
                                                else: self.SuperH[i][j] = D
                                                self.SuperH[j][i] = np.conjugate(self.SuperH[i][j])
                                        # self.log(f"nns: {self.atoms[n_r]} and {self.atoms[n_c]}")

                                # if self.atoms[n_r][0] == "Sn":
                                #     # This is a Hca type block
                                #     if self.atoms[n_c][0] == "Te":  # Hg is the cation
                                #         # print(f"n_r is: {n_r} and s_r is {s_r} and s_c is {s_c}")
                                #         for o_r in self.O:
                                #             for o_c in self.O:
                                #                 # Here we assume that since we are considering only nearest neighbours we have the other type of atoms in teh coulmns
                                #                 i, j, D = self.get_r_c_d(s_r, s_c , n_r, n_c, o_r, o_c)
                                #                 if s_c == s_r: self.SuperH[i][j] = self.Hca[o_r][o_c] + D
                                #                 else: self.SuperH[i][j] = D
                                #         # self.log(f"nns: {self.atoms[n_r]} and {self.atoms[n_c]}")
                                #     if self.atoms[n_c][0] == "H":
                                #         # print(f"n_r is: {n_r} and s_r is {s_r} and s_c is {s_c}")
                                #         for o_r in self.O:
                                #             for o_c in self.H:
                                #                 i, j, D = self.get_r_c_d(s_r, s_c , n_r, n_c, o_r, o_c)
                                #                 if s_c == s_r: self.SuperH[i][j] = self.Hch[o_r][o_c] + D
                                #                 else: self.SuperH[i][j] = D
                                #         # self.log(f"nns: {self.atoms[n_r]} and {self.atoms[n_c]}")

                                if self.atoms[n_r][0] == "H":
                                    # This is a Hha Hc type block
                                    if self.atoms[n_c][0] == "Sn":  # Hg is the cation
                                        # print(f"n_r is: {n_r} and s_r is {s_r} and s_c is {s_c}")
                                        for o_r in self.H:
                                            for o_c in self.O:
                                                # Here we assume that since we are considering only nearest neighbours we have the other type of atoms in teh coulmns
                                                i, j, D = self.get_r_c_d(s_r, s_c , n_r, n_c, o_r, o_c)
                                                # print(i,j)
                                                if s_c == s_r: self.SuperH[i][j] = self.Hha[o_r][o_c] + D
                                                else: self.SuperH[i][j] = D
                                        # self.log(f"nns: {self.atoms[n_r]} and {self.atoms[n_c]}")
                                    # if self.atoms[n_c][0] == "Sn":
                                    #     # print(f"n_r is: {n_r} and s_r is {s_r} and s_c is {s_c}")
                                    #     for o_r in self.H:
                                    #         for o_c in self.O:
                                    #             i, j, D = self.get_r_c_d(s_r, s_c , n_r, n_c, o_r, o_c)
                                    #             if s_c == s_r: self.SuperH[i][j] = self.Hhc[o_r][o_c] + D
                                    #             else: self.SuperH[i][j] = D
                                        # self.log(f"nns: {self.atoms[n_r]} and {self.atoms[n_c]}")
        # self.log("--------- End Nearest neighbours")

        if self.print_hamiltonian:
            self.log(f"\n--------  Printing Full Hamiltonian (Non zero elements) (Before diagonalization):")
            print(f"\n--------  Printing Full Hamiltonian (Before diagonalization):")
            text = ""
            for i in range(len(self.SuperH)):
                # print(f"Row {i}")
                for j in range(len(self.SuperH)):
                    if not self.SuperH[i][j] == 0:
                        text += f"Element {i:<4}{j:<4}: {self.SuperH[i][j]}\n"
            self.log(text)

        self.log(f"\nStarting diagonalization of super matrix......   Time now is: {datetime.now()}")
        print(f"\nStarting diagonalization of super matrix......   Time now is: {datetime.now()}")
        self.eigen_vals, self.eigen_vectors = LA.eig(self.SuperH)
        self.log(f"\nDone.   Time now is: {datetime.now()}")
        print(f"\nDone.   Time now is: {datetime.now()}")

        self.sorted_eigen_vals = self.eigen_vals.copy()
        self.sorted_eigen_vals.sort()
        print(f"\nPrinting values to .out file: {datetime.now()}")

        self.log(f"\n--------  Results: ")
        self.log(f"--------  Eigen values: ")
        for i,v in enumerate(self.sorted_eigen_vals):
            self.log(f"{i+1:<5} = {v:>4}")  # I+1 so that we try to correct for the index


        self.plot()

        self.log(f"\n\nEnd of run >> {datetime.now()}")
        self.end_time = datetime.now()
        self.log(f"Elapased wall time: {self.end_time - self.start_time}")
        self.log(f"Job completed")