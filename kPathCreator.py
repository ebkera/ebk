import numpy as np

# Defined here are some of the high symmetry points in units of 2pi/a
X = [0., 1. ,0., "X"]
L = [1/2, 1/2, 1/2, "L"]
W = [1/2, 1., 0., "W"] 
U = [1/4, 1., 1/4, "U"]
K = [3/4, 3/4, 0., "K"]
G = [0., 0., 0., "$\\Gamma$"]

class kPathCreator():
    def __init__(self):
        self.k_path = []
        self.k_distance = []
        self.highSymPoints_symbols = []
        self.highSymPoints_position = []

    def add_startval(self, start_v):
        """Makes the required K - path to itereate over. Feed in K values without dividing by a"""
        self.highSymPoints_symbols.append(start_v[3])
        self.k_path.append([start_v[0], start_v[1], start_v[2]])
        self.k_distance.append(0)
        self.highSymPoints_position.append(self.k_distance[-1])

    def add_kpath(self, end_v, n_points):
        """Makes the required K - path to itereate over. Feed in K values without dividing by a"""
        # try:
        start_v = self.k_path[-1]
        x = np.linspace(start_v[0],end_v[0],n_points+1)
        y = np.linspace(start_v[1],end_v[1],n_points+1)
        z = np.linspace(start_v[2],end_v[2],n_points+1)
        self.highSymPoints_symbols.append(end_v[3])
        s = self.k_distance[-1]  # getting the current k distance
        for i in range(1,len(x)):
            s = np.sqrt((x[i]-x[i-1])**2 + (y[i]-y[i-1])**2 + (z[i]-z[i-1])**2) + s
            self.k_path.append([x[i], y[i], z[i]])
            self.k_distance.append(s)
        self.highSymPoints_position.append(self.k_distance[-1])
        # except:
        #     print("Start point probably not set!")

    def out_kpath_QE(self, name="kpath"):
        """Prints the kpath in the required format"""
        f = open(f"{name}.kpath","w+")
        f.write(f"#  This file contains the k_path in Quantum Espresso compatible format.\n")
        f.write(f"#  This path contains {len(self.k_distance)} k points.\n")
        f.write(f"#  The path for High Symmetry points  : {self.highSymPoints_symbols}\n")
        f.write(f"#  Locations for High Symmetry points : {self.highSymPoints_position}\n\n")

        for x in range(0,len(self.k_path)):
            f.write(f"{self.k_path[x][0]:.2f}  {self.k_path[x][1]:.2f}  {self.k_path[x][2]:.2f}  {self.k_distance[x]}\n")
        f.close()

if __name__ == "__main__":
    path = kPathCreator()
    path.add_startval(G)
    path.add_kpath(X, 20)
    path.add_kpath(W, 20)
    path.add_kpath(L, 20)
    path.add_kpath(G, 20)
    path.add_kpath(K, 20)
    path.add_kpath(L, 20)
    path.out_kpath_QE()