from ebk.SIESTA.SIESTAOutFileReader import SiestaReadOut
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as cons

class SIESTAcalculatePOT(SiestaReadOut):
    def __init__(self, out_file_name, **kwargs):
        super().__init__(out_file_name)

        self.eps_1 = kwargs.get("eps_1",cons.epsilon_0)
        self.eps_2 = kwargs.get("eps_2",cons.epsilon_0)

        # self.load_quadrupole_moments("/home/ebk/OneDrive_UIC/CQD Research/Analysis/Quadrupole_fitting/TPDAc_F_B3.electrostatics_edit.out")
        self.load_quadrupole_moments()
        # print(self.get_pot_difference())

    def calculate_pointcharge_coordinates(self):
        """This method will calculate the coordinates for the disccrete charge densities"""
        self.point_charge_coordinates = []
        self.point_charge_charge = [1,-1,-1,1]
        self.Q_n = [self.Q_electronic_non_traceless[0][0], self.Q_electronic_non_traceless[1][1], self.Q_electronic_non_traceless[2][2]]
        self.Q_p = [self.Q_ionic_non_traceless[0][0], self.Q_ionic_non_traceless[1][1], self.Q_ionic_non_traceless[2][2]]
        self.q_n = float(self.total_normalized_electronic_charge/2)
        self.q_p = float(self.total_ionic_charge/2)

        def get_pos(Q,q):
            return np.sqrt(Q/q)
        
        positive_point_charge_coordinates = []
        negative_point_charge_coordinates = []
        for i in range(0,3):
            positive_point_charge_coordinates.append(get_pos(self.Q_p[i],self.q_p*2))
            negative_point_charge_coordinates.append(get_pos(self.Q_n[i],self.q_n*2))

        # self.point_charge_coordinates.append([ positive_point_charge_coordinates[0],  positive_point_charge_coordinates[1], -positive_point_charge_coordinates[2]])  #i=1
        # self.point_charge_coordinates.append([ negative_point_charge_coordinates[0],  negative_point_charge_coordinates[1], -negative_point_charge_coordinates[2]])  #i=2
        self.point_charge_coordinates.append([ negative_point_charge_coordinates[0],  negative_point_charge_coordinates[1],  negative_point_charge_coordinates[2]])  #i=3
        self.point_charge_coordinates.append([ positive_point_charge_coordinates[0],  positive_point_charge_coordinates[1],  positive_point_charge_coordinates[2]])  #i=4

        for i,p in enumerate(self.point_charge_coordinates):
            if self.point_charge_charge[i] < 0: 
                Qxz = self.Q_electronic_non_traceless[0][2]
                Qyz = self.Q_electronic_non_traceless[1][2]
            else:
                Qxz = self.Q_ionic_non_traceless[0][2]
                Qyz = self.Q_ionic_non_traceless[1][2]
        
            # Figuring out the sign on the coordinates
            if (Qxz/self.point_charge_charge[i])<0 :
                self.point_charge_coordinates[i][0]=self.point_charge_coordinates[i][0]*(-1)
            if (Qyz/self.point_charge_charge[i])<0:
                self.point_charge_coordinates[i][1]=self.point_charge_coordinates[i][1]*(-1)

        # Adding the other charges onto the front of the list to complete the 4 charges
        self.point_charge_coordinates.insert(0,[-self.point_charge_coordinates[-2][0],-self.point_charge_coordinates[-2][1],-self.point_charge_coordinates[-2][2]])
        self.point_charge_coordinates.insert(0,[-self.point_charge_coordinates[-1][0],-self.point_charge_coordinates[-1][1],-self.point_charge_coordinates[-1][2]])

        # fig,ax = plt.subplots()
        # fig.suptitle(self.SystemLabel)
        # for i,x in enumerate(self.point_charge_coordinates):
        #     if i in [0,3]:s="+"
        #     else : s="."
        #     plt.subplot(2, 1, 1)
        #     plt.plot(x[2],x[0],s,label=f"{i+1}")
        #     plt.xlabel("z")
        #     plt.ylabel("x")

        #     plt.subplot(2, 2, 3)
        #     plt.plot(x[0],x[1],s,label=f"{i+1}")
        #     plt.xlabel("x")
        #     plt.ylabel("y")
        #     plt.subplot(2, 2, 4)
            
        #     plt.plot(x[2],x[1],s,label=f"{i+1}")
        #     plt.xlabel("z")
        #     plt.ylabel("y")
        # plt.legend()
        # plt.show()

    def calculate_coc_quadrupole(self):
        self.coc_quadrupole = np.zeros((3,3))
        self.coc_positive_quadrupole = np.zeros((3,3))
        self.coc_negative_quadrupole = np.zeros((3,3))
        for i,charge in enumerate(self.point_charge_coordinates):
            for row in range(0,3):
                for col in range(0,3):
                    if i == 0 or i == 3:
                        self.coc_positive_quadrupole[row,col] += self.q_p*((charge[row])*(charge[col]))
                    else:
                        self.coc_negative_quadrupole[row,col] += self.q_n*((charge[row])*(charge[col]))
        self.total_coc_quadrupole = self.coc_negative_quadrupole+self.coc_positive_quadrupole

    def get_potential_at_point_from_plate():
        pass

    def get_pot_difference(self):
        a,b = [self.get_initial_cell_vectors()[0][0]*10**(-10),self.get_initial_cell_vectors()[1][1]*10**(-10)]
        c = self.get_initial_cell_vectors()[2][2]
        f = 0.2081943*10**(-10)*10**(-10)  # factor fr converting the quadurupole units to 
        first_term = cons.elementary_charge*f*self.Q_non_traceless[2][2]*(1/self.eps_1 + 1/self.eps_2)/(2*a*b*(np.sqrt(f*self.Q_electronic_non_traceless[2][2]/self.total_normalized_electronic_charge) + np.sqrt(f*self.Q_ionic_non_traceless[2][2]/self.total_ionic_charge)))

        def getA(k1_max,k2_max):
            A = 0
            for k1 in range(1,k1_max):
                for k2 in range(1,k2_max):
                    print("athule")
                    Mk = np.sqrt(k1**2/a**2 + k2**2/b**2)
                    print(np.e)
                    # A+=np.e**

        getA(2,2)

        total = first_term
        # print(a,b)
        return total



