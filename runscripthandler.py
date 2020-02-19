"""
This file creates:
SIRSTAinput files
QE input files
bash scripts for running on CARBON.
bash scripts for running on local machines
This can be used to run multiple jobs for example
As of now this code can only do runs with similar job schedular parameters
"""

import os
from ase.build import bulk
from ase import Atoms
import ase.io
from ebk.QE import QErunfilecreator  # So that we can see how we did it last time
from ebk.SIESTA import SIESTARunFileCreator  # So that we can see how we did it last time

pseudo_database_path = {"cluster":"/usr/local/share/espresso/pseudo",
                        "carbon":"/mnt/c/Users/Eranjan/Desktop/PseudopotentialDatabase",
                        "siva_labs_wsl":"/mnt/c/Users/Eranjan/Desktop/PseudopotentialDatabase",
                        "home_wsl":"/mnt/c/Users/Eranjan/Desktop/PseudopotentialDatabase"
                        }

class RunScriptHandler:
    """All handling of files and script for creating and executing runs is the functionality of this class"""
    def __init__(self, *args, **kwargs):
        """
        kwargs:
            "calc" (string): scf, relax, bands
            "KE_cut" (list): The kinetic energy cutoff
            "identifier" (string): Description of run will be used in file names
            "job_handler" (string): ("slurm", "torque") This is required and will print the right job files for "slurm" or "torque" job handles
            "a0" (list): The lattice constant
            "k" (list of lists whith length 3): The k grid
            "pseudopotentials" (string):
        """
        self.d = f"^"  # Here you can set the desired delimiter
        self.equals = f"+"  # Here you can set the desired symbol for value assigner

        # Gettings args here
            # No args to get currently

        # Setting kwargs here
        # Base Run inits
        self.identifier = kwargs.get("identifier", "run")
        self.job_handler = kwargs.get("job_handler", "torque")
        self.a0 = kwargs.get("a0", [6.6, 6.7, 6.8, 6.9])
        self.KE_cut = kwargs.get("KE_cut", [20, 40, 60, 80, 100])
        self.k = kwargs.get("k", [2])
        self.pseudopotentials = kwargs.get("pseudopotentials", {'Sn': 'Sn_ONCV_PBE_FR-1.1.upf'})
        # self.R = kwargs.get("R", [300])

        # Quantum espresso inits
        self.ntasks = kwargs.get("ntasks", 1)
        self.calc = kwargs.get("calc", "scf")
        self.lspinorb        = True,
        self.noncolin        = True,
        # self.ecutrho         = KE_cut_i*4,
        self.occupations     = 'smearing',
        self.smearing        = 'gaussian',
        self.degauss         = 0.01,
        self.mixing_beta     = 0.7

        # Here goes the job init stuff
        self.walltime_days = 0
        self.walltime_mins = 0
        self.walltime_hours = 2
        self.walltime_secs = 0
        self.nodes = 2
        self.procs = 8

        # Other Initializations
        self.structure = kwargs.get("structure", 1)
        default = Atoms('Sn2', [(0, 0, 0), (0.25, 0.25, 0.25)],  pbc=True)
        self.atoms_object = kwargs.get("atoms_object", default)
        self.all_runs_list = []

    def write_QE_inputfile(self, run_name, KE_cut_i, a0_i, k_i):
        """
        This method creates Quantum espresso input files
        """
        ase.io.write(f"{self.identifier}.in", self.atoms_object, format = "espresso-in", 
                        label           = f"{run_name}",
                        pseudopotentials= self.pseudopotentials,
                        pseudo_dir      = "/mnt/c/Users/Eranjan/Desktop/PseudopotentialDatabase",
                        kpts            = (k_i, k_i, k_i),
                        ecutwfc         = KE_cut_i,
                        calculation     = f"{calc}",
                        lspinorb        = self.lspinorb,
                        noncolin        = self.noncolin,
                        # ecutrho         = KE_cut_i*4,
                        occupations     = self.occupations,
                        smearing        = self.smearing,
                        degauss         = self.degauss,
                        mixing_beta     = self.mixing_beta)

    def write_SIESTA_inputfile(self, run_name, KE_cut_i, a0_i, k_i):
        """
        This method creates SIESTA input files
        """
        pass
        # ase.io.write(f"{self.identifier}.in", self.atoms_object, format = "espresso-in", 
        #                 label           = f"{run_name}",
        #                 pseudopotentials= self.pseudopotentials,
        #                 pseudo_dir      = "/mnt/c/Users/Eranjan/Desktop/PseudopotentialDatabase",
        #                 kpts            = (k_i, k_i, k_i),
        #                 ecutwfc         = KE_cut_i,
        #                 calculation     = f"{calc}",
        #                 lspinorb        = self.lspinorb,
        #                 noncolin        = self.noncolin,
        #                 # ecutrho         = KE_cut_i*4,
        #                 occupations     = self.occupations,
        #                 smearing        = self.smearing,
        #                 degauss         = self.degauss,
        #                 mixing_beta     = self.mixing_beta)

    def get_number_of_calculations(self):
        return (len(self.KE_cut)*len(self.a0)*len(self.k))

    def create_torque_job(self, run_name):
        with open (f"{run_name}.job", "w") as file:
            file.write(f"#!/bin/bash\n")
            file.write(f"#\n")
            file.write(f"#  Basics: Number of nodes, processors per node (ppn), and walltime (hhh:mm:ss)\n")
            file.write(f"#PBS -l nodes={self.nodes}:ppn={self.procs}\n")
            file.write(f"#PBS -l walltime=0:{self.walltime_mins}:00\n")
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
            file.write(f"mpirun -machinefile  $PBS_NODEFILE -np $PBS_NP pw.x -in {run_name}.in -out {run_name}.out\n")

    def create_slurm_job(self, run_name):
            with open (f"{file_name}.job", "w") as file:
                file.write(f"#!/bin/bash\n")
                file.write(f"#SBATCH --job-name={file_name}\n")
                file.write(f"#SBATCH --partition={self.partition}\n")
                file.write(f"#SBATCH --time={self.walltime_days}-{self.walltime_hours}:{self.walltime_mins}:{self.walltime_secs}\n")
                file.write(f"#SBATCH --nodes={self.nodes}\n")
                file.write(f"#SBATCH --ntasks={self.ntasks}\n")
                file.write(f"mpirun -np {self.ntasks} pw.x -in {file_name}.in > {file_name}.out\n")            

    def make_runs(self):
        """This is more Doc strings"""
        # Here the name version for the pseudopotentials is created since we have to have a form that can go on the folder names
        self.PP = ""
        for key,val in self.pseudopotentials:
            self.PP = f"{self.PP}-{val}"
        print(f"PP; {self.PP}")
        if self.structure == None:
            print(f"make_runs: Warning! No structure set. Cannot create run.")
        else:
            for KE_cut_i in self.KE_cut:
                for a0_i in self.a0:
                    for k_i in self.k:
                        # for R_i in self.R:
                        run_name = f"{run_name}KE{self.equals}{KE_cut_i}{self.d}K{self.equals}{k_i}{self.d}R{self.equals}{R_i}{self.d}a{self.equals}{a0_i}{self.d}PP{self.equals}{PP}{self.d}type{self.equals}{calc}"
                        if self.structure == 1:
                            b = a0_i/2.0
                            self.atoms_object.set_cell([(0, b, b), (b, 0, b), (b, b, 0)], scale_atoms=True)
                        else:
                            print("make_runs: Warning! Cannot set cell. Structrue not supported")
                        self.write_QE_inputfile(run_name, KE_cut_i, a0_i, k_i)
                        os.mkdir(f"{run_name}")
                        os.rename(f"{run_name}.in", f"./{run_name}/{run_name}.in")
                        os.rename(f"{run_name}.job", f"./{run_name}/{run_name}.job")
                        self.all_runs_list.append(run_name)

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
        self.pseudopotentials = kwargs.get("pseudopotentials", None)
        self.calc = kwargs.get("calc", None)
        # self.R = kwargs.get("R", None)

        # Here goes the input file parameters stuff
        self.PP = "Sn_ONCV_PBE_FR-1.1.upf"
        self.structure = None

        # Initializations
        self.E_val = []
        self.E_f_val = []
        self.k_val = []
        self.KE_val = []
        self.a0_val = []

    def read_outfiles(self):
        """
        This method reads the out files from the requried directories
        """
        pass
