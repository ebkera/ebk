""" This is is for testing the quadrupole"""
from matplotlib.pyplot import close
import numpy as np

class discrete_quadrupole():
    def __init__(self, CO_dist=1.15908, dipole_from_siesta=[-0.000000, -0.000000, 0.149545]) -> None:
        # Setting the default values here
        self.point_charges = []  # will contain all the charges in this format: [[x,y,z], 1/-1] ([position],charge)

    def set_position(self, position, charge) -> None:
        self.point_charges.append([position,charge])

    def flush_positions(self):
        self.point_charges = []

    def flip_positions_z(self):
        to_append = []
        for point_charge in self.point_charges:
            to_append.append([[point_charge[0][0],point_charge[0][1],-point_charge[0][2]], point_charge[1]])
        self.point_charges.extend(to_append)

    def calculate_quadrupoles(self):
        self.Q = np.zeros((3,3))
        self.Q_non_traceless = np.zeros((3,3))
        self.Q_electronic = np.zeros((3,3))
        self.Q_electronic_non_traceless = np.zeros((3,3))
        self.Q_ionic = np.zeros((3,3))
        self.Q_ionic_non_traceless = np.zeros((3,3))
        self.ionic_charge = 0
        self.electronic_charge = 0
        for point_charge in self.point_charges:
            # print(point_charge)
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
                        self.electronic_charge+= point_charge[1]

                    if point_charge[1] >0:
                        self.Q_ionic[row,col]               += point_charge[1]*(3*(point_charge[0][row])*(point_charge[0][col]) - r2*f)
                        self.Q_ionic_non_traceless[row,col] += point_charge[1]*((point_charge[0][row])*(point_charge[0][col]))
                        self.ionic_charge+= point_charge[1]

    def calculate_quadrupoles_using_quadrupolar_coc(self) -> float:
        coc_electronic = np.sqrt(self.Q_electronic[2][2]/self.electronic_charge)
        coc_ionic = np.sqrt(self.Q_ionic[2][2]/self.ionic_charge)
        d = coc_ionic - coc_electronic
        r = coc_electronic
        Qzz = (self.ionic_charge/2)*d*(2*r+d)
        return Qzz

    def convert_to_debye(self):
        unit_factor_Debye = 1/0.2081943 #(same as c.elementary_charge*10**(-10)/(3.33564*10**(-30)))
        # print("printing the Q and debye")
        # print(self.Q)
        # print(self.Q_non_traceless)
        # print(self.Q_debye)
        self.Q = self.Q*unit_factor_Debye
        self.Q_non_traceless = self.Q_non_traceless*unit_factor_Debye
        self.Q_electronic = self.Q_electronic*unit_factor_Debye
        self.Q_electronic_non_traceless = self.Q_electronic_non_traceless*unit_factor_Debye
        self.Q_ionic = self.Q_ionic*unit_factor_Debye
        self.Q_ionic_non_traceless = self.Q_ionic_non_traceless*unit_factor_Debye

    def plot_charge_profile(self, x_range = 15, title_text="Z axis", save_to="discrete_point_charges_vs_coc"):
        import matplotlib.pyplot as plt 
        print("Plotting charge profile. Please wait...")
        resolution= 500
        x = list(np.linspace(-x_range, x_range, resolution))
        y_n = np.linspace(0, 0, resolution)
        y_p = np.linspace(0, 0, resolution)
        for p_i, p in enumerate(self.point_charges):
            for x_i in x:
                if p[0][2] <= x_i:
                    if p[1]>0:y_p[x.index(x_i)] += p[1]
                    if p[1]<0:y_n[x.index(x_i)] += p[1]
                    break
        plt.plot(x,y_n, label = "negative part")
        plt.plot(x,y_p, label = "positive part")
        plt.xlabel("The z axis $\AA$")
        plt.ylabel("Charge in electrons")
        # plt.title(f"Charge density profile {title_text}")
        plt.legend()
        plt.savefig(f"{save_to}/{title_text}.pdf")
        # plt.show()
        plt.close()
        # fig(close)


