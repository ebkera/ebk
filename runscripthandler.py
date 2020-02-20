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
from ebk.calculation_set import Calculation_set
import shutil


class RunScriptHandler():
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
            "atoms_object" (atoms object): This should be without setting the cell since that will be done with every iteration
            "structure" (int): vlaue will determine the cell
                1: fcc structure with a0 as the lattice constant

        """
        self.d = f"^"  # Here you can set the desired delimiter
        self.equals = f"="  # Here you can set the desired symbol for value assigner

        # Gettings args here
            # No args to get currently

        # Setting kwargs here
        # Base Run inits
        self.identifier       = kwargs.get("identifier", "run")
        self.job_handler      = kwargs.get("job_handler", "torque")
        self.a0               = kwargs.get("a0", [6.6, 6.7, 6.8, 6.9])
        self.KE_cut           = kwargs.get("KE_cut", [20, 40, 60, 80, 100])
        self.k                = kwargs.get("k", [2])
        test = {'Sn':'Sn_ONCV_PBE_FR-1.1.upf','B':'B_ONCV_PBE_FR-1.1.upf'}
        self.pseudopotentials = kwargs.get("pseudopotentials", test)
        self.pseudo_dir       = kwargs.get("pseudo_dir", None)
        self.calculator       = kwargs.get("calculator", "espresso")
        # self.R = kwargs.get("R", [300])

        # Quantum espresso inits
        self.ntasks = kwargs.get("ntasks", 20)
        self.calc = kwargs.get("calc", "scf")
        self.lspinorb        = kwargs.get("lspinorbit", False)
        self.noncolin        = kwargs.get("noncolin", False)
        # self.ecutrho         = kwargs.get("KE_cut_i*4,
        self.occupations     = kwargs.get("occupations",'smearing')
        self.smearing        = kwargs.get("smearing",'gaussian')
        self.degauss         = kwargs.get("degauss", 0.01)
        self.mixing_beta     = kwargs.get("mixing_beta", 0.7)
        self.xc              = kwargs.get("xc", "pbe")
        self.structure_type  = kwargs.get("structure_type", "bulk")
        self.Title           = kwargs.get("Title",'EthaneDithiol')
        self.prefix          = kwargs.get("prefix",'E2D')
        self.restart_mode    = kwargs.get("restart_mode",'from_scratch')
        self.disk_io         = kwargs.get("disk_io",'default')
        self.verbosity       = kwargs.get("verbosity",'high')
        self.lkpoint_dir     = kwargs.get("lkpoint_dir", False)
        self.etot_conv_thr   = kwargs.get("etot_conv_thr", 1.0e-6)
        self.forc_conv_thr   = kwargs.get("forc_conv_thr", 1.0e-4)
        self.outdir          = kwargs.get("outdir", './')
        self.pseudo_dir      = kwargs.get("pseudo_dir", False)

        # Here goes the job init stuff
        self.walltime_days   = kwargs.get("walltime_days", 2)
        self.walltime_mins   = kwargs.get("walltime_mins", 0)
        self.walltime_hours  = kwargs.get("walltime_hours", 2)
        self.walltime_secs   = kwargs.get("walltime_secs", 0)
        self.nodes           = kwargs.get("nodes", 2)
        self.procs           = kwargs.get("procs", 20)
        self.partition       = kwargs.get("partition", "cluster")

        # Other Initializations
        self.structure = kwargs.get("structure", 1)
        default = Atoms('Sn2', [(0, 0, 0), (0.25, 0.25, 0.25)],  pbc=True)
        self.atoms_object = kwargs.get("atoms_object", default)
        self.all_runs_list = []

    def set_pseudo_dir(self, machine):
        """
        Sets the pseudo_dir according to the machine
        """
        pseudo_database_path = {"cluster":"/usr/local/share/espresso/pseudo",
                        "carbon":"/mnt/c/Users/Eranjan/Desktop/PseudopotentialDatabase",
                        "siva_labs_wsl":"/mnt/c/Users/Eranjan/Desktop/PseudopotentialDatabase",
                        "home_wsl":"/mnt/c/Users/Eranjan/Desktop/PseudopotentialDatabase"
                        }
        self.pseudo_dir = pseudo_database_path[machine]

    def write_QE_inputfile(self, run_name, KE_cut_i, a0_i, k_i):
        """
        This method creates Quantum espresso input files
        """
        ase.io.write(f"{self.identifier}.in", self.atoms_object, format = "espresso-in", 
                        label           = f"{run_name}",
                        pseudopotentials= self.pseudopotentials,
                        # if self.pseudo_dir == None
                        pseudo_dir      = self.pseudo_dir,
                        kpts            = (k_i, k_i, k_i),
                        ecutwfc         = KE_cut_i,
                        calculation     = f"{self.calc}",
                        lspinorb        = self.lspinorb,
                        noncolin        = self.noncolin,
                        # ecutrho         = KE_cut_i*4,
                        occupations     = self.occupations,
                        smearing        = self.smearing,
                        degauss         = self.degauss,
                        mixing_beta     = self.mixing_beta,
                        Title           = self.Title,
                        prefix          = self.prefix,
                        restart_mode    = self.restart_mode,
                        disk_io         = self.disk_io,
                        verbosity       = self.verbosity,
                        lkpoint_dir     = self.lkpoint_dir,
                        etot_conv_thr   = self.etot_conv_thr,
                        forc_conv_thr   = self.forc_conv_thr,
                        outdir          = self.outdir,
                        )

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
        with open (f"{self.identifier}.job", "w") as file:
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
        os.rename(f"{self.identifier}.job", f"./{run_name}/{self.identifier}.job")

    def create_slurm_job(self, run_name):
        with open (f"{self.identifier}.job", "w") as file:
            file.write(f"#!/bin/bash\n")
            file.write(f"#SBATCH --job-name={run_name}\n")
            file.write(f"#SBATCH --partition={self.partition}\n")
            file.write(f"#SBATCH --time={self.walltime_days}-{self.walltime_hours}:{self.walltime_mins}:{self.walltime_secs}\n")
            file.write(f"#SBATCH --nodes={self.nodes}\n")
            file.write(f"#SBATCH --ntasks={self.procs}\n")
            # file.write(f"#SBATCH --mail-user=erathnayake@sivananthanlabs.us\n")
            # file.write(f"#SBATCH --mail-type=ALL\n")
            file.write(f"mpirun -np {self.ntasks} pw.x < {self.identifier}.in > {self.identifier}.out\n")
        os.rename(f"{self.identifier}.job", f"./{run_name}/{self.identifier}.job")

    def make_runs(self):
        """This is more Doc strings"""
        # Here the name version for the pseudopotentials is created since we have to have a form that can go on the folder names
        self.PP = ""
        self.specie = ""
        for key, val in self.pseudopotentials.items():
            self.PP = f"{val}-{self.PP}"
            self.specie = f"{key}-{self.specie}"

        if self.structure == None:
            print(f"make_runs: Warning! No structure set. Cannot create run.")
        else:
            for KE_cut_i in self.KE_cut:
                for a0_i in self.a0:
                    for k_i in self.k:
                        # for R_i in self.R:  # This has been disables for now
                        R_i = KE_cut_i*4
                        run_name = f"{self.identifier}{self.d}Calc{self.equals}{self.calculator}{self.d}Struct{self.equals}{self.structure_type}{self.d}Specie{self.equals}{self.specie}{self.d}KE{self.equals}{KE_cut_i}{self.d}K{self.equals}{k_i}{self.d}R{self.equals}{R_i}{self.d}a{self.equals}{a0_i}{self.d}PP{self.equals}{self.PP}{self.d}type{self.equals}{self.calc}"
                        if self.structure == 0:
                            # cell has been set from outside
                            pass
                        elif self.structure == 1:
                            # An fcc cell that scales with the lattice constant
                            b = a0_i/2.0
                            self.atoms_object.set_cell([(0, b, b), (b, 0, b), (b, b, 0)], scale_atoms=True)
                        else:
                            print("make_runs: Warning! Cannot set cell. Structrue not supported")
                        self.write_QE_inputfile(run_name, KE_cut_i, a0_i, k_i)
                        if os.path.exists(run_name):
                            shutil.rmtree(run_name)
                            print("make_runs: Warning! Path exists!! Overwriting")
                        os.mkdir(f"{run_name}")
                        os.rename(f"{self.identifier}.in", f"./{run_name}/{self.identifier}.in")
                        self.all_runs_list.append(run_name)
                        
                        # Creating jobs
                        if self.job_handler == "torque":
                            self.create_torque_job(run_name)
                        elif self.job_handler == "slurm":
                            self.create_slurm_job(run_name)
                        else:
                            print(f"make_runs: Unrecognized job_handler! Job files not created")

    def create_bash_file(self):
        """
        This script creates bash files so that you can run a batch of the runs that need to be done
        """
        print(f"create_bash_file: job_handler is set to: {self.job_handler}")
        bash_file = open("run.sh", "w+")
        bash_file.write(f"#!/bin/bash\n\n")
        bash_file.write(f"dir_list=(")
        for x in self.all_runs_list:
            bash_file.write(f" {x}")
        bash_file.write(")\n")
        bash_file.write('for dir in "${dir_list[@]}"\n')
        bash_file.write(f"do\n")
        bash_file.write(f'  cd $dir\n')
        bash_file.write(f'  dos2unix *job\n')

        if self.job_handler == "torque":
            bash_file.write(f'  qsub *job\n')
        if self.job_handler == "slurm":
            bash_file.write(f'  sbatch *job\n')
        bash_file.write(f'  cd ..\n')
        bash_file.write(f'done\n')
        bash_file.close()


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

if __name__ == "__main__":
    """This is used as an example as to how we use this file."""
    pass