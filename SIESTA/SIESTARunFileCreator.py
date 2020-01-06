import os
import shutil
import sys
import subprocess
from ebk.kPathCreator import *
from ebk import A2Bohr

class SIESTARunCreator:
    def __init__(self, system_name):
        """
        This function initializes the object. Default initializations are
        |self.job_manager [options: "torque", "slurm"]
        """
        # Initialization for job
        self.wall_time = [0, 30, 0]
        self.nodes = 4
        self.procs = 8
        self.job_manager = "torque" #options: (torque, slurm)
        # Initialization for run
        self.system_name = system_name
        self.k_set = [30,32,34]
        self.ke_set = [180, 200, 220]
        self.r_set = [300, 400]
        self.a_set = [6.7, 6.8, 6.9]
        self.k_points_number = 0
        self.excorr = "pbe"     # "pbe"   = "sla+pw+pbx+pbc"    = PBE, "pz"    = "sla+pz"            = Perdew-Zunger LDA
        self.PP_name = "Sn.UPF"
        # self.celldm = [0, 12.394714, 0, 0, 0, 0, 0] # celldm variables. index 0 is not used and index 1,6 corresponds to celldm1,6
        self.lspinorb = False
        self.supercell_file = "default"
        self.atomic_species = [["Sn", "118.71", "Sn.UPF"]]
        self.lat_const = [6.5]
        self.celldm = [0, 12.394714, 0, 0, 0, 0, 0] # celldm variables. index 0 is not used and index 1,6 corresponds to celldm1,6
        self.cell_parameters = False
        self.set_of_files = []
        self.ibrav = 2
        self.nat = 2
        self.ntyp = 1
        self.dirname = "Run"

    def make_name(self):
        for k in self.k_set:
            for ke in self.ke_set:
                for r in self.r_set:
                    for a in self.a_set:
                        name = f"{self.system_name}_SIESTA_K{k}_KE{ke}_R{r}_a{a:.2f}"
        return name

if __name__ == "__main__":
    Sn_run = QERunCreator("Sn")
    path = kPathCreator()
    path.add_startval(G)
    path.add_kpath(X, 20)
    path.add_kpath(W, 20)
    path.add_kpath(L, 20)
    path.add_kpath(G, 20)
    path.add_kpath(K, 20)
    path.add_kpath(L, 20)
    path.out_kpath_QE()
    Sn_run.k_points_number = len(path.k_distance)
    Sn_run.excorr = 'auto'
    # Sn_run.excorr = 'pz'
    Sn_run.k_set = [10]
    Sn_run.ke_set = [40]
    Sn_run.r_set = [300]
    Sn_run.create_run(30, together = True)