import os
import shutil
import sys
import subprocess
from ebk.kPathCreator import *
from ebk import A2Bohr

class QERunCreator:
    def __init__(self, system_name, k_set = [30,32,34], ke_set = [180, 200, 220], r_set = [300, 400]):
        """
        This function initializes the object but usually 
        """
        self.system_name = system_name
        self.k_points_number = 0
        self.k_set = k_set
        self.ke_set = ke_set
        self.r_set = r_set
        self.excorr = "pbe"     # "pbe"   = "sla+pw+pbx+pbc"    = PBE, "pz"    = "sla+pz"            = Perdew-Zunger LDA
        self.PP_name = "Sn.UPF"
        self.nodes = 4
        self.procs = 8
        # self.celldm = [0, 12.394714, 0, 0, 0, 0, 0] # celldm variables. index 0 is not used and index 1,6 corresponds to celldm1,6
        self.lspinorb = False
        self.supercell_file = "default"
        self.atomic_species = [["Sn", "118.71", "Sn.UPF"]]
        self.lat_const = [6.5]
        self.celldm = [0, 12.394714, 0, 0, 0, 0, 0] # celldm variables. index 0 is not used and index 1,6 corresponds to celldm1,6
        self.set_of_files = []

    def make_name(self, k, ke, r, bands, a = False):
        name = f"{self.system_name}_QE_K{k}_KE{ke}_R{r}"
        if a != False:
            name = f"{name}_a{a:.2f}"

        # Saving the set of file names generated for later use
        self.set_of_files.append(name)
        return name

    def jobCreator(self, k, ke, r, walltime_mins, bands, dirname, a):
        name = self.make_name(k, ke, r, bands, a)
        if bands == True:
            file_name = f"{name}.bands"
        else:
            file_name = f'{name}.scf'
        with open (f"{file_name}.job", "w") as file:
            file.write(f"#!/bin/bash\n")
            file.write(f"#\n")
            file.write(f"#  Basics: Number of nodes, processors per node (ppn), and walltime (hhh:mm:ss)\n")
            file.write(f"#PBS -l nodes={self.nodes}:ppn={self.procs}\n")
            file.write(f"#PBS -l walltime=0:{walltime_mins}:00\n")
            file.write(f"#PBS -N {file_name}\n")
            file.write(f"#PBS -A cnm66441\n")
            file.write(f"#\n")
            file.write(f"#  File names for stdout and stderr.  If not set here, the defaults\n")
            file.write(f"#  are <JOBNAME>.o<JOBNUM> and <JOBNAME>.e<JOBNUM>\n")
            file.write(f"#PBS -o job.out\n")
            file.write(f"#PBS -e job.err\n")
            file.write("\n")
            file.write(f"# Send mail at begin, end, abort, or never (b, e, a, n). Default is 'a'.\n")
            file.write(f"#PBS -m bea erathnayake@sivananthanlabs.us\n")
            file.write("\n")
            file.write(f"# change into the directory where qsub will be executed\n")
            file.write(f"cd $PBS_O_WORKDIR\n")
            file.write("\n")
            file.write(f"# start MPI job over default interconnect; count allocated cores on the fly.\n")
            file.write(f"mpirun -machinefile  $PBS_NODEFILE -np $PBS_NP pw.x -in {file_name}.in > {file_name}.out\n")
        shutil.move(f"{file_name}.job", f"./{dirname}/{file_name}.job")

    def k_file_reader(self):
        with open ("kpath.kpath",'r') as k_path_file:
            self.k_path = k_path_file.read()

    def infileCreator(self, k, ke, r, bands, dirname, a):
        # We need to have the below lines since the prefix has to not have the .scf / .bands parts
        name = self.make_name(k, ke, r, bands, a)
        if bands == True:
            file_name = f"{name}.bands"
        else:
            file_name = f'{name}.scf'
        with open (f"{file_name}.in", "w") as file:
            file.write(f"&control\n")
            if bands == False:
                file.write(f"    calculation     = 'scf'\n")
            else:
                file.write(f"    calculation     = 'bands'\n")
            file.write(f"    verbosity       = 'high'\n")
            if self.excorr != "auto":
                file.write(f"    input_dft       = '{self.excorr}'\n")
            file.write(f"    prefix          = '{name}'\n")
            file.write(f"    wf_collect      = .false.\n")
            file.write(f"    pseudo_dir      = './'\n")
            file.write(f"    outdir          = './'\n")
            file.write(f"/\n")
            file.write(f"&system\n")
            file.write(f"    ibrav           = 2\n")
            if a != False:
                # Works only if a lattice parameter has been set.
                file.write(f"    celldm(1)       = {a}\n")
            else:
                for x in range(1,6):
                    if self.celldm[x] != 0:
                        file.write(f"    celldm({x})       = {self.celldm[x]}\n")
            file.write(f"    nat             = 2\n")
            file.write(f"    ntyp            = 1\n")
            file.write(f"    ecutwfc         = {ke}\n")
            if bands == False:
                pass
            else:
                file.write(f"    nbnd            = 30\n")
            if self.lspinorb == True:
                file.write(f"    lspinorb        = .true.\n")
                file.write(f"    noncolin        = .true.\n")
            file.write(f"    ecutrho         = {r}\n")
            file.write(f"    occupations     = 'smearing'\n")
            file.write(f"    smearing        = 'gaussian'\n")
            file.write(f"    degauss         = 0.01\n")
            file.write(f"/\n")
            file.write(f"&electrons\n")
            file.write(f"    mixing_beta     = 0.7\n")
            file.write(f"/\n")
            file.write(f"ATOMIC_SPECIES\n")
            for specie in self.atomic_species:
                file.write(f"{specie[0]}   {specie[1]}   {specie[2]}\n")
            file.write(f"ATOMIC_POSITIONS alat \n")
            # Here we have to make srue that the new requrired supercell (basis) is loaded into the file
            if self.supercell_file == "default":
                file.write(f"Sn        0.000000000   0.000000000   0.000000000\n")
                file.write(f"Sn        0.250000000   0.250000000   0.250000000\n")
            else:
                try:
                    with open (f"{self.supercell_file}", "r") as file2:
                        data = file2.read()
                        data = data.split("\n")
                        for line in data:
                            # x = line.split()
                            # print
                            if len(line.split()) == 4:
                                file.write(f"{line}\n")
                except:
                    print(f"QERunCreator.inFileCreator: Could not open file {self.supercell_file}. Writing default structure (diamond) with prefix Sn")
                    file.write(f"Sn        0.000000000   0.000000000   0.000000000\n")
                    file.write(f"Sn        0.250000000   0.250000000   0.250000000\n")

            if bands == False:
                file.write(f"K_POINTS automatic\n{k} {k} {k} 1 1 1\n")
            else:
                file.write(f"\n")
                file.write(f"K_POINTS tpiba\n#  tpiba = k-points in units of 2pi/a\n#  number of k-points\n")
                file.write(f"{self.k_points_number}\n")
                file.write(f"#  kx, ky, kz, wk  (wk=symmetry weight is ignored in bands calculations)\n")
                file.write(self.k_path)
        shutil.move(f"{file_name}.in", f"./{dirname}/{file_name}.in")

    def create_run(self, walltime, together):
        self.k_file_reader()
        if together == True:
            dirname = "Run"
        else:
            dirname = f"{self.system_name}_QE_K{k}_KE{ke}_R{r}"
        for k in self.k_set:
            for ke in self.ke_set:
                for r in self.r_set:
                    try:
                        os.mkdir(dirname)
                    except:
                        pass
                    for a in self.lat_const:
                        b = A2Bohr(a)
                        self.jobCreator(k, ke, r, walltime, True, dirname, b)
                        self.jobCreator(k, ke, r, walltime, False, dirname, b)
                        self.infileCreator(k, ke, r, False, dirname, b)
                        self.infileCreator(k, ke, r, True, dirname, b)

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
