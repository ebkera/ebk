""" This is is for testing the quadrupole"""

from email.policy import default
from turtle import position
from pandas import array
""" This is the file that makes the game of the game """
import numpy as np

class OCCO_quadrupole():
    def __init__(self) -> None:
        # Setting the default values here
        self.CO_dist = 1.15908
            # THis is in angstroms
        self.CO_dipole_from_siesta = np.array([-0.000000, -0.000000, 0.149545])
            # This dipole value is unrelaxed which is because relaxation would make it loose is CO relaxation value
            # The values are in Debye
        self.CO_dipole_in_eA = 0.2081943*self.CO_dipole_from_siesta
        self.charge_on_a_pointcharge = self.CO_dipole_in_eA[2]/self.CO_dist
            # In number of electrons
        self.point_charges = []  # will contain all the charges in this format: [[x,y,z], 1/-1] ([position],charge)

    def set_positions_for_zz(self, CC_distance:float) -> None:
        self.point_charges.append([[0,0,-CC_distance/2],1*self.charge_on_a_pointcharge])
        self.point_charges.append([[0,0,CC_distance/2],1*self.charge_on_a_pointcharge])
        self.point_charges.append([[0,0,-CC_distance/2-self.CO_dist],-1*self.charge_on_a_pointcharge])
        self.point_charges.append([[0,0,CC_distance/2+self.CO_dist],-1*self.charge_on_a_pointcharge])
        
    def calculate_quadrupoles(self) -> list:
        self.Q = np.zeros((3,3))
        self.Q_non_traceless = np.zeros((3,3))
        self.Q_electronic = np.zeros((3,3))
        self.Q_electronic_non_traceless = np.zeros((3,3))
        self.Q_ionic = np.zeros((3,3))
        self.Q_ionic_non_traceless = np.zeros((3,3))
        for point_charge in self.point_charges:
            r_vec = point_charge[0]
            r2 = np.dot(r_vec, r_vec)
            for col in range(0,3):
                for row in range(0,3):
                    if col == row: f=1
                    else: f=0
                    self.Q[row,col]               += point_charge[1]*(3*(point_charge[0][row])*(point_charge[0][col]) - r2*f)
                    self.Q_non_traceless[row,col] += point_charge[1]*((point_charge[0][row])*(point_charge[0][col]))

                    if point_charge[1] <0:
                        self.Q_electronic[row,col]               += point_charge[1]*(3*(point_charge[0][row])*(point_charge[0][col]) - r2*f)
                        self.Q_electronic_non_traceless[row,col] += point_charge[1]*((point_charge[0][row])*(point_charge[0][col]))

                    if point_charge[1] >0:
                        self.Q_ionic[row,col]               += point_charge[1]*(3*(point_charge[0][row])*(point_charge[0][col]) - r2*f)
                        self.Q_ionic_non_traceless[row,col] += point_charge[1]*((point_charge[0][row])*(point_charge[0][col]))

    def convert_to_debye(self):
        unit_factor_Debye = 1/0.2081943 #(same as c.elementary_charge*10**(-10)/(3.33564*10**(-30)))
        self.Q_debye = self.Q*unit_factor_Debye
        self.Q_non_traceless = self.Q_non_traceless*unit_factor_Debye
        self.Q_electronic_debye = self.Q_electronic*unit_factor_Debye
        self.Q_electronic_non_traceless = self.Q_electronic_non_traceless*unit_factor_Debye
        self.Q_ionic_debye = self.Q_ionic*unit_factor_Debye
        self.Q_ionic_non_traceless = self.Q_ionic_non_traceless*unit_factor_Debye

    def print_values_to_console(self):
        # Some testing prints
        print("\nCalculated values:\n")
        print("dipoles from siesta (DA, eA): ", self.CO_dipole_from_siesta, self.CO_dipole_in_eA)
        unit_factor_Debye = 1/0.2081943 #(same as c.elementary_charge*10**(-10)/(3.33564*10**(-30)))
        print(f"Q in eA          : \n{self.Q}")
        print(f"Qtraceless in eA : \n{self.Q_non_traceless}")
        Q = self.Q*unit_factor_Debye
        Q_non_traceless = self.Q_non_traceless*unit_factor_Debye
        print(f"Q in D.A         : \n{Q}")
        print(f"Qtraceless in D.A: \n{Q_non_traceless}")

