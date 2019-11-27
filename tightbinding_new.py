import numpy as np
from numpy import linalg as LA
import math
# import pandas as pd
import matplotlib.pyplot as plt
# from ebk.kPathCreator import kPathCreator
from ebk.BandPlotter import BandPlotter
from ebk.kPathCreator import *

# def convert_eV_to_J(E_in):
#     return E_in * 1.60218e-19
    
# def convert_J_to_eV(E_in):
#     return E_in * 6.242e+18

# def calculate_lmn(vector):
#     """ This function caluclates lmn values for each R vector"""
#     # I don't use this function as of now
#     R = np.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)
#     l = vector[0]/R
#     m = vector[1]/R
#     n = vector[2]/R
#     # print (f"calculate_lmn: {l,m,n}")  # For testing purposes
#     return l,m,n

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

class TightBindingCalculation():
    def __init__(self, name):
        # Mapping of numbers and orbitals, sites, and spins
        self.i_s = 0
        self.i_px = 1
        self.i_py = 2
        self.i_pz = 3
        self.i_dxy = 4
        self.i_dyz = 5
        self.i_dzx = 6
        self.i_dz = 7  #dx2-y2
        self.i_d3 = 8  #d3z^2-r^2
        self.i_S = 9  #s star
        self.i_anion = 0
        self.i_cation = 10
        self.i_up = 0
        self.i_down = 20

        # Some initializations
        self.name = name
        self.a = 1 * 6.453e-10
        self.k_path = []
        self.k_distance = []
        self.eigen_vals = []
        self.eigen_vectors = []
        self.H = np.zeros((20,20),dtype=complex)
        self.cations_as_neighbours = np.array([[0.25,0.25,0.25],[0.25,-0.25,-0.25],[-0.25,0.25,-0.25],[-0.25,-0.25,0.25]])*self.a
        # self.anions_as_neighbours = np.array([[-0.25,-0.25,-0.25],[-0.25,0.25,0.25],[0.25,-0.25,0.25],[0.25,0.25,-0.25]])*self.a
        self.so = False

    def load_parameters(self):
        """This method will eventually be able to load from file. Since we dont need it now we can just do a quick assignment of the parameters"""
        # We save all the parameters here into a dictionary
        self.parameter = {}

        if self.so == True:
            self.so_a = 0.375000  # This is in eV
            self.so_c = 0.465000  # This is in eV
        elif self.so == False:
            self.so_a = 0.0  # This is in eV
            self.so_c = 0.0  # This is in eV

        # self.parameter.update({"so_a": self.so_a})
        # self.parameter.update({"so_c": self.so_c})

        # # Defining tight bindging paramters here _ca only for diagonal elements (remember we are using the ac block matrix as the base and then getting ca just by complex conjugate transpose)
        # # The anion energies
        # self.parameter.update({"Es_a": -5.9819})
        # self.parameter.update({"Ep_a": 3.5820})
        # self.parameter.update({"Edxy_a": 13.1023})
        # self.parameter.update({"Edx2my2_a": 13.1023})
        # self.parameter.update({"ES_a": 19.4220})

        # # The cation energies
        # self.parameter.update({"Es_c": -0.4028})
        # self.parameter.update({"Ep_c": 6.3853})
        # self.parameter.update({"Edxy_c": 13.1023})
        # self.parameter.update({"Edx2my2_c": 13.1023})
        # self.parameter.update({"ES_c": 19.4220})

        # # Now we load the other parameters
        # # The diagonal parameters
        # self.parameter.update({"Vsssigma": -1.6187})
        # self.parameter.update({"Vddsigma": -0.529629})
        # self.parameter.update({"Vppsigma": 3.166827})
        # self.parameter.update({"Vpppi": -0.945694})
        # self.parameter.update({"Vddpi": 2.424709})
        # self.parameter.update({"Vdddelta": -1.064207})
        # self.parameter.update({"VSSsigma": -3.6761})

        # # The anion to cation off-diagonal parameters
        # self.parameter.update({"Vspsigma_ac": 1.085069})
        # self.parameter.update({"Vsdsigma_ac": -0.525896})
        # self.parameter.update({"Vpdsigma_ac": -1.789915})
        # self.parameter.update({"Vpdpi_ac": 1.406422})
        # self.parameter.update({"VsS_ac": 1.9927})
        # self.parameter.update({"VSpsigma_ac": 1.175059})
        # self.parameter.update({"VSdsigma_ac": 0.485896})

        # # The cation to anion off-diagonal parameters
        # self.parameter.update({"Vspsigma_ca": 2.014492})
        # self.parameter.update({"Vsdsigma_ca": -1.067102})
        # self.parameter.update({"Vpdsigma_ca": -0.653612})
        # self.parameter.update({"Vpdpi_ca": 1.657517})
        # self.parameter.update({"VsS_ca": -0.242580})
        # self.parameter.update({"VSpsigma_ca": 1.405375})
        # self.parameter.update({"VSdsigma_ca": 0.696627})

        # This is for HgTe
        self.parameter.update({"so_a": 2*self.so_a})
        self.parameter.update({"so_c": self.so_c})

        # Defining tight bindging paramters here _ca only for diagonal elements (remember we are using the ac block matrix as the base and then getting ca just by complex conjugate transpose)
        # The anion energies
        self.parameter.update({"Es_a": -10.040161})
        self.parameter.update({"Ep_a": 1.580003})
        self.parameter.update({"Edxy_a": 10.139959})
        self.parameter.update({"Edx2my2_a": 13.145395})
        self.parameter.update({"ES_a": 12.611213})

        # The cation energies
        self.parameter.update({"Es_c": -1.502103})
        self.parameter.update({"Ep_c": 5.929255})
        self.parameter.update({"Edxy_c": 15.108978})
        self.parameter.update({"Edx2my2_c": 15.431086})
        self.parameter.update({"ES_c": 14.801158})

        # Now we load the other parameters
        # The diagonal parameters
        self.parameter.update({"Vsssigma": -0.904384})
        self.parameter.update({"Vddsigma": -0.529629})
        self.parameter.update({"Vppsigma": 3.166827})
        self.parameter.update({"Vpppi": -0.945694})
        self.parameter.update({"Vddpi": 2.424709})
        self.parameter.update({"Vdddelta": -1.064207})
        self.parameter.update({"VSSsigma": -1.570513})

        # The anion to cation off-diagonal parameters
        self.parameter.update({"Vspsigma_ac": 1.085069})
        self.parameter.update({"Vsdsigma_ac": -0.525896})
        self.parameter.update({"Vpdsigma_ac": -1.789915})
        self.parameter.update({"Vpdpi_ac": 1.406422})
        self.parameter.update({"VsS_ac": 0.357261})
        self.parameter.update({"VSpsigma_ac": 1.175059})
        self.parameter.update({"VSdsigma_ac": 0.485896})

        # The cation to anion off-diagonal parameters
        self.parameter.update({"Vspsigma_ca": 2.014492})
        self.parameter.update({"Vsdsigma_ca": -1.067102})
        self.parameter.update({"Vpdsigma_ca": -0.653612})
        self.parameter.update({"Vpdpi_ca": 1.657517})
        self.parameter.update({"VsS_ca": -0.242580})
        self.parameter.update({"VSpsigma_ca": 1.405375})
        self.parameter.update({"VSdsigma_ca": 0.696627})

        # Save new list with all paramters converted to J
        # self.parameterJ = {key: value * 1.60218e-19 for (key,value) in self.parameter.items()}
        self.parameterJ = {key: value for (key,value) in self.parameter.items()}

    def Vsp(self, direction):
        if direction == "ac":
            return self.parameterJ["Vspsigma_ac"]/np.sqrt(3)
        elif direction == "ca":
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

    def add_kpath(self,kpath):
        """Makes the required K - path to itereate over. Feed in K values without dividing by a"""
        self.k_path = []
        for point in kpath:
            self.k_path.append([x*(2*np.pi)/self.a for x in point])

    def print_kpath(self):
        f = open(f"{self.name}.k","w+")
        for x in self.k_path:
            f.write(f"{x[0]}\t{x[1]}\t{x[2]}\n")
        f.close()

    def create_Hamiltonian(self):
        """indeces 0-9 are anions and 10-19 are cations."""
        big_block = [0,10]  # anions and cations
        small_block = list(range(0,10))
        self.e_file= open(f"{self.name}.elements","w+")
        eig_file= open(f"{self.name}.eig","w+")
        check = open(f"{self.name}.check","w+")
        check2 = open(f"{self.name}.check2","w+")

        for k_point in self.k_path:
            self.e_file.write(f"Kpoint {k_point}###############################################\n")
            eig_file.write(f"Kpoint {k_point}###############################################\n")
            check.write(f"Kpoint {k_point}###############################################\n")
            check2.write(f"Kpoint {k_point}###############################################\n")

            for x1 in big_block:
                for x2 in big_block:
                    for y1 in small_block:
                        for y2 in small_block:
                            row = x1 + y1
                            col = x2 + y2
                            element = self.function_chooser(row, col, k_point)
                            # Writing to elements file
                            self.e_file.write(f"row:{row}, column: {col} val:{element}\n")
                            if row > 9 and col < row-10:
                                self.e_file.write(f"Orig: {self.H[col][row]} Now(conjugated): {np.conjugate(element)} Diff: {self.H[col][row]-np.conjugate(element)}\n")
                            self.H[row][col] = element

            # Testing the matrix for errors
            for row2 in range(0,len(self.H)):
                for col2 in range(0,len(self.H)):
                    if self.H[col2][row2]-np.conjugate(self.H[row2][col2]) != 0:

                        check2.write(f"row:{row2}, col:{col2}\n")
                        check2.write(f"Orig: {self.H[col2][row2]} Now(conjugated): {np.conjugate(self.H[row2][col2])} Diff: {self.H[col2][row2]-np.conjugate(self.H[row2][col2])}\n")

            if self.so == False:
                # Diagonalize and calculate the eigen values for no Spin orbit coupling
                w,v = LA.eig(self.H)

            elif self.so == True:
                # Diagonalize and calculate the eigen values for no Spin orbit coupling
                self.SuperH = np.zeros((40,40),dtype=complex)

                # Lets first create the super matrix made up of two submatrices of self.H. 
                # The off diagonal submatrices are zero and the diagonal are the calculated self.H matrices
                for row in range(0,len(self.H)):
                    for column in range(0,len(self.H)):
                        self.SuperH[row][column] = self.H[row][column]
                        self.SuperH[row+20][column+20] = self.H[row][column]

                # Now we look at the spin orbit coupling elements
                # Lets do the up-up sub-matrix first
                self.SuperH[self.i_up+self.i_anion+self.i_py][self.i_up+self.i_anion+self.i_px] += self.parameterJ["so_a"]*1j
                self.SuperH[self.i_up+self.i_anion+self.i_px][self.i_up+self.i_anion+self.i_py] += -self.parameterJ["so_a"]*1j
                self.SuperH[self.i_up+self.i_cation+self.i_py][self.i_up+self.i_cation+self.i_px] += self.parameterJ["so_c"]*1j
                self.SuperH[self.i_up+self.i_cation+self.i_px][self.i_up+self.i_cation+self.i_py] += -self.parameterJ["so_c"]*1j
                # Lets do the down-down sub-matrix next
                self.SuperH[self.i_down+self.i_anion+self.i_py][self.i_down+self.i_anion+self.i_px] += -self.parameterJ["so_a"]*1j
                self.SuperH[self.i_down+self.i_anion+self.i_px][self.i_down+self.i_anion+self.i_py] += self.parameterJ["so_a"]*1j
                self.SuperH[self.i_down+self.i_cation+self.i_py][self.i_down+self.i_cation+self.i_px] += -self.parameterJ["so_c"]*1j
                self.SuperH[self.i_down+self.i_cation+self.i_px][self.i_down+self.i_cation+self.i_py] += self.parameterJ["so_c"]*1j
                # Lets do the up-down sub-matrix next
                self.SuperH[self.i_up+self.i_anion+self.i_pz][self.i_down+self.i_anion+self.i_px] += -self.parameterJ["so_a"]
                self.SuperH[self.i_up+self.i_anion+self.i_px][self.i_down+self.i_anion+self.i_pz] += self.parameterJ["so_a"]
                self.SuperH[self.i_up+self.i_anion+self.i_pz][self.i_down+self.i_anion+self.i_py] += self.parameterJ["so_a"]*1j
                self.SuperH[self.i_up+self.i_anion+self.i_py][self.i_down+self.i_anion+self.i_pz] += -self.parameterJ["so_a"]*1j
                self.SuperH[self.i_up+self.i_cation+self.i_pz][self.i_down+self.i_cation+self.i_px] += -self.parameterJ["so_c"]
                self.SuperH[self.i_up+self.i_cation+self.i_px][self.i_down+self.i_cation+self.i_pz] += self.parameterJ["so_c"]
                self.SuperH[self.i_up+self.i_cation+self.i_pz][self.i_down+self.i_cation+self.i_py] += self.parameterJ["so_c"]*1j
                self.SuperH[self.i_up+self.i_cation+self.i_py][self.i_down+self.i_cation+self.i_pz] += -self.parameterJ["so_c"]*1j
                # Lets do the down-up sub-matrix next
                self.SuperH[self.i_down+self.i_anion+self.i_px][self.i_up+self.i_anion+self.i_pz] += -self.parameterJ["so_a"]
                self.SuperH[self.i_down+self.i_anion+self.i_pz][self.i_up+self.i_anion+self.i_px] += self.parameterJ["so_a"]
                self.SuperH[self.i_down+self.i_anion+self.i_py][self.i_up+self.i_anion+self.i_pz] += -self.parameterJ["so_a"]*1j
                self.SuperH[self.i_down+self.i_anion+self.i_pz][self.i_up+self.i_anion+self.i_py] += self.parameterJ["so_a"]*1j
                self.SuperH[self.i_down+self.i_cation+self.i_px][self.i_up+self.i_cation+self.i_pz] += -self.parameterJ["so_c"]
                self.SuperH[self.i_down+self.i_cation+self.i_pz][self.i_up+self.i_cation+self.i_px] += self.parameterJ["so_c"]
                self.SuperH[self.i_down+self.i_cation+self.i_py][self.i_up+self.i_cation+self.i_pz] += -self.parameterJ["so_c"]*1j
                self.SuperH[self.i_down+self.i_cation+self.i_pz][self.i_up+self.i_cation+self.i_py] += self.parameterJ["so_c"]*1j
                w,v = LA.eig(self.SuperH)
            # w = convert_J_to_eV(w)

            # for row1 in range(0,len(self.SuperH)):
            #     for col1 in range(0,len(self.SuperH)):
            #         if self.SuperH[col1][row1]-np.conjugate(self.SuperH[row1][col1]) != 0:
            #             check.write(f"row:{row1}, col:{col1}\n")
            #             check.write(f"Orig: {self.SuperH[col1][row1]} Now(conjugated): {np.conjugate(self.SuperH[row1][col1])} Diff: {self.SuperH[col1][row1]-np.conjugate(self.SuperH[row1][col1])}\n")

            w.sort()
            eig_file.write(f"{str(w)}\n")
            self.eigen_vals.append(w)
            self.eigen_vectors.append(v)

            # Writing to Hamil file
            H_file= open(f"{self.name}.hamil","w+")
            for x in self.H:
                for y in x:
                    H_file.write(f'{y:.2e}\t')
                H_file.write("\n")
            H_file.close()
        self.e_file.close()
        eig_file.close()
        check.close()
        check2.close()


    def function_chooser(self, row, col, k_point):
        if row == col:
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
                else:
                    print("function_chooser: Fatal Error - Element not found")

        elif row <= 9 and col > 9:
            # These are the anion to cation elements
            col = col - 10  # So we can still refer to the rows and columns as from 0-9 to retain our sanity
            # next let's select the diagonal and upper triangle
            neighbours_ac = self.cations_as_neighbours
            # neighbours_ac = self.anions_as_neighbours
            # neighbours_ca = self.anions_as_neighbours
            # neighbours_ca = self.cations_as_neighbours

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
            # These are the cations to anion elements
            # They will be the conjugate transpose of the anion to cation elements
            # Therefore we have to make sure that this will by done after the row<=9 and col > 9 set is done with
            return np.conjugate(self.H[col][row])
        else:
            return 0. + 0.j

if __name__ == "__main__":
    path = kPathCreator()
    path.add_startval(L)
    path.add_kpath(G, 20)
    path.add_kpath(X, 20)
    path.add_kpath(W, 20)
    path.add_kpath(G, 20)
    path.add_kpath(K, 20)
    # path.out_kpath_QE()

    HgTe = TightBindingCalculation("HgTe")
    HgTe.so = True
    HgTe.load_parameters()
    HgTe.add_kpath(path.k_path)

    HgTe.create_Hamiltonian()
    HgTe.print_kpath()
    HgTe.eigen_vals = np.transpose(HgTe.eigen_vals)
    plotter= BandPlotter(path.k_distance, HgTe.eigen_vals, path.highSymPoints_position, path.highSymPoints_symbols)
    plotter.title = "TightBinding HgTe"
    plotter.file_name = "BandStructure"
    plotter.vlines = True
    plotter.set_y_range = True
    plotter.ylim_high = 15
    plotter.ylim_low = -15
    plotter.same_band_colour = False
    plotter.plot()