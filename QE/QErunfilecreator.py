import os
import shutil
import sys
import subprocess
from ebk.kPathCreator import *

class QERunCreator:
    def __init__(self, system_name):
        self.system_name = system_name
        
    def make_name(self, k, ke, r, bands):
        return f"{self.system_name}_QE_K{k}_KE{ke}_R{r}"

    def jobCreator(self, k, ke, r, bands, dirname):
        name = self.make_name(k, ke, r, bands)
        file_name = name
        if bands == True:
            file_name = f"{name}.bands"
        with open (f"{file_name}.job", "w") as file:
            file.write(f"#!/bin/bash\n")
            file.write(f"#\n")
            file.write(f"#  Basics: Number of nodes, processors per node (ppn), and walltime (hhh:mm:ss)\n#PBS -l nodes=5:ppn=8\n")
            file.write(f"#PBS -l walltime=0:10:00\n")
            if bands == True:
                file.write(f"#PBS -N {name}.bands\n")
            else:
                file.write(f"#PBS -N {name}\n")
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
            if bands == True:
                file.write(f"mpirun -machinefile  $PBS_NODEFILE -np $PBS_NP pw.x -in {name}.bands.in > {name}.bands.out\n")
            else:
                file.write(f"mpirun -machinefile  $PBS_NODEFILE -np $PBS_NP pw.x -in {name}.in > {name}.out\n")
        shutil.move(f"{file_name}.job", f"./{dirname}/{file_name}.job")

    def k_file_reader(self):
        with open ("kpath.kpath",'r') as k_path_file:
            self.k_path = k_path_file.read()

    def infileCreator(self, k, ke, r, bands, dirname):
        name = self.make_name(k, ke, r, bands)
        file_name = name
        if bands == True:
            file_name = f"{name}.bands"
        with open (f"{file_name}.in", "w") as file:
            file.write(f"&control\n")
            if bands == False:
                file.write(f"    calculation     = 'scf'\n")
            else:
                file.write(f"    calculation     = 'bands'\n")
            file.write(f"    verbosity       = 'high'\n")
            file.write(f"    prefix          = '{name}'\n")
            file.write(f"    wf_collect      = .false.\n")
            file.write(f"    pseudo_dir      = './'\n")
            file.write(f"    outdir          = './'\n")
            file.write(f"/\n")
            file.write(f"&system\n")
            file.write(f"    ibrav           = 2\n")
            file.write(f"    celldm(1)       = 12.394714\n")
            file.write(f"    nat             = 2\n")
            file.write(f"    ntyp            = 1\n")
            file.write(f"    ecutwfc         = {ke}\n")
            if bands == False:
                pass
            else:
                file.write(f"    nbnd            = 30\n")
            file.write(f"    ecutrho         = {r}\n")
            file.write(f"    occupations     = 'smearing'\n")
            file.write(f"    smearing        = 'gaussian'\n")
            file.write(f"    degauss         = 0.01\n")
            file.write(f"/\n")
            file.write(f"&electrons\n")
            file.write(f"    mixing_beta     = 0.7\n")
            file.write(f"/\n")
            file.write(f"ATOMIC_SPECIES\n")
            file.write(f" Sn 118.71  Sn.UPF\n")
            file.write(f"ATOMIC_POSITIONS alat \n")
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


    def runcreator(self, k_set, ke_set, r_set, together):
        self.k_file_reader()
        if together == True:
            dirname = "Run"
        else:
            dirname = f"{self.system_name}_QE_K{k}_KE{ke}_R{r}"
        for k in k_set:
            for ke in ke_set:
                for r in r_set:
                    try:
                        os.mkdir(dirname)
                    except:
                        pass
                    self.jobCreator(k, ke, r, True, dirname)
                    self.jobCreator(k, ke, r, False, dirname)
                    self.infileCreator(k, ke, r, False, dirname)
                    self.infileCreator(k, ke, r, True, dirname)

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
    k_set = [30,31,32]
    ke_set = [200, 250, 300, 320]
    r_set = [300, 400]
    Sn_run.runcreator(k_set, ke_set, r_set, together = True)