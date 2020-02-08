"""This file creates bash scripts for running on CARBON. This can be used to run multiple jobs for example"""

import os
from ase.build import bulk
from ase import Atoms
import ase.io

class Runscriptcreator:
    """Doc string goes here"""
    def __init__(*args, **kwargs):
        """Doc string goes here"""
        # Here goes the PBS init stuff
        self.walltime_mins = 30
        self.nodes = 2
        self.procs = 8
        # Here goes the other stuff
        self.PP = "Sn_ONCV_PBE_FR-1.1.upf"
        self.pseudopotentials = {'Sn': 'Sn_ONCV_PBE_FR-1.1.upf'}
        self.a0 = kwargs.get("a0", [6.6, 6.7, 6.8, 6.9])
        self.KE_cut = kwargs.get("KE_cut", [20, 40, 60, 80, 100])
        self.E = []
        self.k = kwargs.get("k", [2])
        # self.R = kwargs.get("R", [300])
        self.calc = f"scf"
        d = f"^"  # Here you can set the desired delimiter
        equals = f"="
        self.structure = None

    def get_number_of_calculations(self):
        return (self.KE_cut.len()*self.a0.len().self.k.len()*self.R.len())

    def create(self):
        """This is more Doc strings"""
        if self.structure == None:
            print(f"No structure set")
        else:
            for KE_cut_i in self.KE_cut:
                for a0_i in self.a0:
                    for k_i in self.k:
                        for R_i in self.R:
                            run_name = f"QE{d}KE{KE_cut_i}{d}K{k_i}{d}R{R_i}{d}a{a0_i}{d}PP={PP}{d}calc={calc}"
                            self.structure.set_cell([(b, 0, 0), (0, b, 0), (0, 0, b)], scale_atoms=True)
                            ase.io.write(f"{run_name}.in", bulk, format = "espresso-in", 
                                            label           = f"{run_name}",
                                            pseudopotentials= pseudopotentials,
                                            pseudo_dir      = "/mnt/c/Users/Eranjan/Desktop/PseudopotentialDatabase",
                                            kpts            = (k_i, k_i, k_i),
                                            ecutwfc         = KE_cut_i,
                                            calculation     = f"{calc}",
                                            lspinorb        = True,
                                            noncolin        = True,
                                            # ecutrho         = KE_cut_i*4,
                                            occupations     = 'smearing',
                                            smearing        = 'gaussian',
                                            degauss         = 0.01,
                                            mixing_beta     = 0.7)

                            with open (f"{run_name}.job", "w") as file:
                                file.write(f"#!/bin/bash\n")
                                file.write(f"#\n")
                                file.write(f"#  Basics: Number of nodes, processors per node (ppn), and walltime (hhh:mm:ss)\n")
                                file.write(f"#PBS -l nodes={nodes}:ppn={procs}\n")
                                file.write(f"#PBS -l walltime=0:{walltime_mins}:00\n")
                                file.write(f"#PBS -N {run_name}\n")
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
                                file.write(f"mpirun -machinefile  $PBS_NODEFILE -np $PBS_NP pw.x -in {run_name}.in > {run_name}.out\n")
                            os.mkdir(f"{run_name}")
                            os.rename(f"{run_name}.in", f"./{run_name}/{run_name}.in")
                            os.rename(f"{run_name}.job", f"./{run_name}/{run_name}.job")