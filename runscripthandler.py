"""
This file creates:
SIRSTAinput files
QE input files
bash scripts for running on CARBON.
bash scripts for running on local machines
This can be used to run multiple jobs for example
"""

import os
from ase.build import bulk
from ase import Atoms
import ase.io

pseudo_database_path = {"cluster":"/usr/local/share/espresso/pseudo",
                        "carbon":"/mnt/c/Users/Eranjan/Desktop/PseudopotentialDatabase",
                        "siva_labs_wsl":"/mnt/c/Users/Eranjan/Desktop/PseudopotentialDatabase",
                        "home_wsl":"/mnt/c/Users/Eranjan/Desktop/PseudopotentialDatabase",
                        }

class RunScriptHandler:
    """Doc string goes here"""
    def __init__(*args, **kwargs):
        """Doc string goes here"""
        self.d = f"^"  # Here you can set the desired delimiter
        self.equals = f"+"  # Here you can set the desired symbol for value assigner

        # Gettings args here

        # Getting kwargs here
        self.a0 = kwargs.get("a0", [6.6, 6.7, 6.8, 6.9])
        self.KE_cut = kwargs.get("KE_cut", [20, 40, 60, 80, 100])
        self.k = kwargs.get("k", [2])
        self.pseudopotentials = kwargs.get("pseudos", {'Sn': 'Sn_ONCV_PBE_FR-1.1.upf'})
        self.calc = kwargs.get("calc", "scf")


        # Here goes the PBS init stuff
        self.walltime_mins = 30
        self.nodes = 2
        self.procs = 8

        # Here goes the input file parameters stuff
        self.PP = "Sn_ONCV_PBE_FR-1.1.upf"

        # self.R = kwargs.get("R", [300])
        d = f"^"  # Here you can set the desired delimiter
        equals = f"="
        self.structure = None

        # Initializations
        self.E = []


    def get_number_of_calculations(self):
        return (self.KE_cut.len()*self.a0.len().self.k.len()*self.R.len())

    def create_bash():
        pass

    def create(self):
        """This is more Doc strings"""
        if self.structure == None:
            print(f"No structure set")
        else:
            for KE_cut_i in self.KE_cut:
                for a0_i in self.a0:
                    for k_i in self.k:
                        for R_i in self.R:
                            run_name = f"QE{d}KE{KE_cut_i}{self.d}K{k_i}{d}R{R_i}{self.d}a{a0_i}{self.d}PP={self.PP}{self.d}calc={self.calc}"
                            self.atoms_object = Atoms('Sn2', [(0, 0, 0), (0.25, 0.25, 0.25)],  pbc=True)
                            b = a0_i/2.0
                            if self.structure[1].len() == 1:
                                self.atoms_object.set_cell([(0, b, b), (b, 0, b), (b, b, 0)], scale_atoms=True)
                            ase.io.write(f"{run_name}.in", self.atoms_object, format = "espresso-in", 
                                            label           = f"{run_name}",
                                            pseudopotentials= self.pseudopotentials,
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

class Read_outfiles():
    """This method should read all out files of a given type (sesta/qe) and read the vlaues like total energies"""
    def __init__(*args, **kwargs):
        """Doc string goes here"""
        self.d = f"^"  # Here you can set the desired delimiter
        self.equals = ["+", "="]  # Here you can set the desired symbol for value assigner it can also be a list of all possible values

        # Gettings args here

        # Getting kwargs here
        self.a0 = kwargs.get("a0", None)
        self.KE_cut = kwargs.get("KE_cut", None)
        self.k = kwargs.get("k", None)
        self.pseudopotentials = kwargs.get("pseudos", None})

        # Here goes the input file parameters stuff
        self.PP = "Sn_ONCV_PBE_FR-1.1.upf"

        # self.R = kwargs.get("R", [300])
        self.calc = f"scf"
        d = f"^"  # Here you can set the desired delimiter
        equals = f"="
        self.structure = None

        # Initializations
        self.E = []


    def read_outfiles(self):
        """This method reads the out files from the requried directories
        pass