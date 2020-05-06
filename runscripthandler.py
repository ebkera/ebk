"""
This file creates:
SIESTAinput files
QE input files
bash scripts for running on torque.
bash scripts for running on local machines
This can be used to run multiple jobs for example
As of now this code can only do runs with similar job scheduler parameters
"""

import os
from ase.build import bulk
from ase import Atoms
import ase.io
from ebk.QE import QErunfilecreator  # So that we can see how we did it last time
from ebk.SIESTA import SIESTARunFileCreator  # So that we can see how we did it last time
import shutil

# Here stored are all the possible values for these variables. They are common for both RunScriptHandler and Read_outfiles classes
identifier       = []
job_handler      = []
a0               = []
KE_cut           = []
k                = []
pseudopotentials = []
pseudo_dir       = []
calculator       = []
structure_type   = []
xc               = []
calculation      = []

class RunScriptHandler():
    """All handling of files and script for creating and executing runs is the functionality of this class"""
    def __init__(self, *args, **kwargs):
        """
        kwargs:
            "calc" (string): scf, relax, bands
            "KE_cut" (list): The kinetic energy cutoff
            "identifier" (string): Description of run will be used in file names
            "job_handler" (string): ("slurm", "torque", "era_pc", "era_ubuntu", "sl_laptop") This is required and will print the right job files for "slurm" or "torque" job handles
            "a0" (list of floats): The lattice constant
            "k" (list of lists whith length 3 or list of ints): The k grid
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
        self.pseudopotentials = kwargs.get("pseudopotentials", {'Sn':'Sn_ONCV_PBE_FR-1.1.upf'})
        self.pseudo_dir       = kwargs.get("pseudo_dir", False)
        self.calculator       = kwargs.get("calculator", "QE")
        self.structure_type   = kwargs.get("structure_type", "bulk")
        self.xc               = kwargs.get("xc", "pbe")
        self.calculation      = kwargs.get("calculation", "scf")
        self.path             = kwargs.get("path", "GXWLGKL")
        self.density          = kwargs.get("density", 30)
        self.k_path           = {"path":self.path, "density": self.density}
        self.R                = kwargs.get("R", [None])

        # Quantum espresso inits some other inits that need to be only set if explicitly given can be found below this.
        self.espresso_inputs = {"pseudopotentials": self.pseudopotentials,
                                "calculation"     : self.calculation,
                                "lspinorb"        : kwargs.get("lspinorb", False),
                                "noncolin"        : kwargs.get("noncolin", False),
                                "occupations"     : kwargs.get("occupations",'smearing'),
                                "diagonalization" : kwargs.get("diagonalization",'david'),
                                "smearing"        : kwargs.get("smearing",'gaussian'),
                                "degauss"         : kwargs.get("degauss", 0.01),
                                "mixing_beta"     : kwargs.get("mixing_beta", 0.7),
                                "Title"           : kwargs.get("Title",'Sn'),
                                "prefix"          : kwargs.get("prefix",'Sn'),
                                "restart_mode"    : kwargs.get("restart_mode",'from_scratch'),
                                "disk_io"         : kwargs.get("disk_io",'default'),
                                "wf_collect"      : kwargs.get("wf_collect", False),
                                "verbosity"       : kwargs.get("verbosity",'high'),
                                "lkpoint_dir"     : kwargs.get("lkpoint_dir", False),
                                "etot_conv_thr"   : kwargs.get("etot_conv_thr", 1.0e-3),
                                "forc_conv_thr"   : kwargs.get("forc_conv_thr", 1.0e-3),
                                "outdir"          : kwargs.get("outdir", './'),
                                "path"            : self.path,
                                "density"         : self.density,  # This is an ASE command for input files for Quantum Espresso
                                "electron_maxstep": kwargs.get("electron_maxstep", 200)
                                }

        # Here are all initializations of the self.espresso_inputs variable that should be set only if explicitly given by user
        if "nbnd" in kwargs:
            self.espresso_inputs.update({"nbnd"            : kwargs.get("nbnd", 40)})


        if self.pseudo_dir != False:
            # Since if not set we want the value to be the default value Thereby reading in machine defaults and not appearing in the .in file
            # pseduo_dir is initialized in the base run init section
            # Doing it this way without haveing it in the quantum espresso inits makes this not appear in the .in file if not set.
            # Thereby reading in machine defaults
            self.espresso_inputs.update({"pseudo_dir" : self.pseudo_dir})

        # Here goes the job init stuff and recources allocation
        self.walltime_days   = kwargs.get("walltime_days", 2)
        self.walltime_mins   = kwargs.get("walltime_mins", 0)
        self.walltime_hours  = kwargs.get("walltime_hours", 2)
        self.walltime_secs   = kwargs.get("walltime_secs", 0)
        self.nodes           = kwargs.get("nodes", 2)  # the number of cores
        self.procs           = kwargs.get("procs", 8)  # number of processesors per core
        self.partition       = kwargs.get("partition", "cluster")  # The partition that the job will run on
        self.ntasks          = kwargs.get("ntasks", 20)  # number of threads in total
        self.npool           = kwargs.get("npool", 1)  # The number of pools of k points per proc (has to be an integer). This is a Quantum espresso parameter and will only work with QE. 

        # Other Initializations
            # For structure: An fcc cell that scales with the lattice constant = 1
        self.structure = kwargs.get("structure", 1)
        default = Atoms('Sn2', [(0, 0, 0), (0.25, 0.25, 0.25)],  pbc=True)
        self.atoms_object = kwargs.get("atoms_object", default)
        self.all_runs_list = []

        # Default executable paths are set here
        self.executable_path = {
                        "cluster":"",
                        "torque":"",
                        "siva_labs_wsl":"/mnt/c/Users/Eranjan/Desktop/PseudopotentialDatabase",
                        "era_pc":"/mnt/c/Users/Eranjan/Desktop/Quantum_Expresso/qe-6.4.1/bin/",
                        "era_ubuntu":"/home/era/Downloads/qe-6.5/bin/",  # the pw - "/home/era/Downloads/qe-6.5/bin"
                        "sl_laptop": "/mnt/c/Users/erathnayake/Desktop/qe-6.5/bin/"
                        }

    def set_pseudo_dir(self, machine):
        """
        Sets the pseudo_dir according to the machine
        """
        pseudo_database_path = {"cluster":"/usr/local/share/espresso/pseudo",
                        "torque":"../PseudopotentialDatabase",
                        "siva_labs_wsl":"/mnt/c/Users/Eranjan/Desktop/PseudopotentialDatabase",
                        "era_pc":"/mnt/c/Users/Eranjan/Desktop/PseudopotentialDatabase",
                        "era_ubuntu":"../PseudopotentialDatabase",  # the pw - "/home/era/Downloads/qe-6.5/bin"
                        "sl_laptop": "/mnt/c/Users/erathnayake/Desktop/PseudopotentialDatabase"
                        }
        self.espresso_inputs.update({"pseudo_dir" : pseudo_database_path[machine]})



    def set_pseudopotentials(self, pseudos):
        """
        Sets the pseudopotentials
        """
        self.pseudopotentials = pseudos
        self.espresso_inputs.update({"pseudopotentials": pseudos})

    def write_QE_inputfile(self, run_name, KE_cut_i, R_i, a0_i, k_i):
        """
        This method creates Quantum espresso input files
        """
        # Here we update aditional run specific stuff
        self.espresso_inputs.update({"label" : f"{run_name}"})
        self.espresso_inputs.update({"ecutwfc" : KE_cut_i})
        if type(k_i) == list:
            self.espresso_inputs.update({"kpts" : (k_i[0], k_i[1], k_i[2])})
        else:
            self.espresso_inputs.update({"kpts" : (k_i, k_i, k_i)})
        if R_i != None: self.espresso_inputs.update({"ecutrho" : R_i})

        if self.calculation == "bands":
            # First we deal with the scf run
            self.espresso_inputs.update({"calculation" : "scf"})
            ase.io.write(f"{self.identifier}.in", self.atoms_object, format = "espresso-in", **self.espresso_inputs)
            # Then here we take care of the bands run
            self.espresso_inputs.update({"calculation" : "bands"})
            self.espresso_inputs.update({"kpts" : self.k_path})
            ase.io.write(f"{self.identifier}.bands.in", self.atoms_object, format = "espresso-in", **self.espresso_inputs)
            os.rename(f"{self.identifier}.bands.in", f"./{run_name}/{self.identifier}.bands.in")
            os.rename(f"{self.identifier}.in", f"./{run_name}/{self.identifier}.in")
        else:
            ase.io.write(f"{self.identifier}.in", self.atoms_object, format = "espresso-in", **self.espresso_inputs)
            os.rename(f"{self.identifier}.in", f"./{run_name}/{self.identifier}.in")

    def write_SIESTA_inputfile(self, run_name, KE_cut_i, a0_i, k_i):
        """
        This method creates SIESTA input files
        """
        pass

    def get_number_of_calculations(self):
        return (len(self.KE_cut)*len(self.a0)*len(self.k)*len(self.R))

    def create_torque_job(self):
        """
        This is the torque job creator and it specifically has PBS lines that torque can understand. 
        The difference between this and the generic job creator is that the generic job creator does not rely on torque.
        """
        # Creating the scf run for bands runs if self.calculation bands
        if self.calculation == "bands":
            self.calculation = "scf"
            self.create_torque_job()
            self.calculation = "bands"

        with open (f"{self.identifier}.{self.calculation}.job", "w+") as file_torque:
            file_torque.write(f"#!/bin/bash\n")
            file_torque.write(f"# Submit jobs from explicitly specified directories;\n")
            file_torque.write(f"# stern, 2020-02-18 - Edited Eranjan\n")
            file_torque.write(f"\n")
            file_torque.write(f'shopt -s extglob	# handle "+()" patterns\n')
            file_torque.write(f"\n")
            file_torque.write(f"# Let's use a shell loop to read a list of tasks from a 'Here-Document' (at the\n")
            file_torque.write(f"# end of the loop).  See also: http://www.tldp.org/LDP/abs/html/here-docs.html\n")
            file_torque.write(f"while read dir\n")
            file_torque.write(f"do\n")
            file_torque.write(f"    # Skip empty lines and comments\n")
            file_torque.write(f'    [[ $dir == +(""|"#"*) ]] && continue\n')
            file_torque.write(f"    # Let's use some basic shell variable string operations to remove the\n")
            file_torque.write(f"    # leading and trailing parts of the dir names that are the same, isolating\n")
            file_torque.write(f"    # the changing middle section as a useful short name.\n")
            file_torque.write("    job_name=${dir#*KE=}\n")
            file_torque.write("    job_name=${job_name%^a*}\n")
            file_torque.write(f"    # Use another here-doc to read the job script, to obviate individual files\n")
            file_torque.write(f"    # in each data directory.\n")
            file_torque.write(f"\n")
            file_torque.write(f"    # A here-doc with leading '-' will get leading TABs removed.  (unwise to\n")
            file_torque.write(f"    # use, however, if your editor munges TABs.)\n")
            file_torque.write(f"\n")
            file_torque.write(f"    # Double-quote $dir to avoid parameter substitution.  (This is not strictly\n")
            file_torque.write(f"    # necessary here, though, because the chars '^+' are not special -- in the\n")
            file_torque.write(f"    # circumstances used here.)\n")
            file_torque.write(f'    qsub -w "$PWD/$dir" -N "$job_name" <<-END_JOB_SCRIPT\n')
            file_torque.write(f"\n")
            file_torque.write(f"    #!/bin/bash\n")
            file_torque.write(f"    #  Basics: Number of nodes, processors per node (ppn), and walltime (hhh:mm:ss)\n")
            file_torque.write(f"    #PBS -l nodes={self.nodes}:ppn={self.procs}\n")
            file_torque.write(f"    #PBS -l walltime={self.walltime_hours}:{self.walltime_mins}:{self.walltime_secs}\n")
            # file.write(f"    #PBS -N {run_name}\n")
            file_torque.write(f"    #PBS -A cnm66441\n")
            file_torque.write(f"    #  File names for stdout and stderr.  If not set here, the defaults\n")
            file_torque.write(f"    #  are <JOBNAME>.o<JOBNUM> and <JOBNAME>.e<JOBNUM>\n")
            # file.write(f"    #PBS -o job.out\n")
            file_torque.write(f"    #PBS -e $PWD/$dir/job.err\n")
            file_torque.write(f"    #PBS -j eo\n")
            file_torque.write(f"    # Send mail at begin, end, abort, or never (b, e, a, n). Default is 'a'.\n")
            file_torque.write(f"    #PBS -m bea\n")
            file_torque.write("\n")
            file_torque.write(f"    # change into the directory where qsub will be executed\n")
            file_torque.write(f"    cd \$PBS_O_WORKDIR\n\n")
            file_torque.write(f"    # use a per-job lineup of modules; stern\n")
            file_torque.write(f"    module purge\n")
            file_torque.write(f"    module load intel\n")
            file_torque.write(f"    module load openmpi/1.10/intel-17\n")
            file_torque.write(f"    module load quantum-espresso/5.4/openmpi-1.10\n")
            file_torque.write(f"    module list\n\n")
            # file.write(f"    # start MPI job over default interconnect; count allocated cores on the fly.\n")
            # file.write(f"    mpirun -machinefile  $PBS_NODEFILE -np $PBS_NP pw.x < {run_name}.in > {run_name}.out\n")
            if self.calculation == "bands":
                file_torque.write(f"    mpirun pw.x < {self.identifier}.bands.in > {self.identifier}.bands.out\n")
            else:
                file_torque.write(f"    mpirun -np {self.ntasks} pw.x -npool {self.npool} < {self.identifier}.in > {self.identifier}.out\n")
            file_torque.write(f"END_JOB_SCRIPT\n")
            file_torque.write(f"\n")
            file_torque.write(f"done <<'END_TASKLIST'\n")
            file_torque.write(f"    # Single quoting the limit string 'EOT' will pass strings without shell variable and execution expansion.\n")
            file_torque.write(f"    # Comments and empty line are fine because we explicitly skip them.\n\n")
            for run in self.all_runs_list:
                file_torque.write(f"    {run}\n")
            file_torque.write(f"END_TASKLIST\n")
        # os.rename(f"{self.identifier}.job", f"./{run_name}/{self.identifier}.job")

    def create_job(self):
        """
        This is the most generic job creator. Will create all jobs
        """
        # Creating the scf run for bands runs if self.calculation bands
        if self.calculation == "bands":
            self.calculation = "scf"
            self.create_job()
            self.calculation = "bands"
        with open (f"{self.identifier}.{self.calculation}.job", "w+") as file_torque:
            file_torque.write(f"#!/bin/bash\n")
            file_torque.write(f"# Submit jobs from explicitly specified directories;\n")
            file_torque.write(f"# stern, 2020-03-18 - Edited Eranjan\n")
            file_torque.write(f"\n")
            file_torque.write(f'shopt -s extglob	# handle "+()" patterns\n')
            file_torque.write(f"\n")
            file_torque.write(f"# Let's use a shell loop to read a list of tasks from a 'Here-Document' (at the\n")
            file_torque.write(f"# end of the loop).  See also: http://www.tldp.org/LDP/abs/html/here-docs.html\n")
            file_torque.write(f"while read dir\n")
            file_torque.write(f"do\n")
            file_torque.write(f"    # Skip empty lines and comments\n")
            file_torque.write(f'    [[ $dir == +(""|"#"*) ]] && continue\n')
            file_torque.write(f"    # Let's use some basic shell variable string operations to remove the\n")
            file_torque.write(f"    # leading and trailing parts of the dir names that are the same, isolating\n")
            file_torque.write(f"    # the changing middle section as a useful short name.\n")
            file_torque.write("    job_name=${dir#*KE=}\n")
            file_torque.write("    job_name=${job_name%^a*}\n")
            file_torque.write(f"    # Use another here-doc to read the job script, to obviate individual files\n")
            file_torque.write(f"    # in each data directory.\n")
            file_torque.write(f"\n")
            file_torque.write(f"    # A here-doc with leading '-' will get leading TABs removed.  (unwise to\n")
            file_torque.write(f"    # use, however, if your editor munges TABs.)\n")
            file_torque.write(f"\n")
            file_torque.write(f"    # Double-quote $dir to avoid parameter substitution.  (This is not strictly\n")
            file_torque.write(f"    # necessary here, though, because the chars '^+' are not special -- in the\n")
            file_torque.write(f"    # circumstances used here.)\n")
            file_torque.write(f'    cd "$PWD/$dir"\n')
            file_torque.write(f"\n")
            file_torque.write(f"    #!/bin/bash\n")
            file_torque.write("\n")
            if self.calculation == "bands":
                file_torque.write(f"    mpirun {self.executable_path[self.job_handler]}pw.x < {self.identifier}.bands.in | tee {self.identifier}.bands.out\n")
            else:
                file_torque.write(f"    mpirun -np {self.ntasks} {self.executable_path[self.job_handler]}pw.x -npool {self.npool} < {self.identifier}.in | tee {self.identifier}.out\n")
            file_torque.write("    cd .. \n")
            file_torque.write(f"\n")
            file_torque.write(f"done <<'END_TASKLIST'\n")
            file_torque.write(f"    # Single quoting the limit string 'EOT' will pass strings without shell variable and execution expansion.\n")
            file_torque.write(f"    # Comments and empty line are fine because we explicitly skip them.\n\n")
            for run in self.all_runs_list:
                file_torque.write(f"    {run}\n")
            file_torque.write(f"END_TASKLIST\n")

    def create_slurm_job(self, run_name):
        with open (f"{self.identifier}.job", "w+") as file:
            file.write(f"#!/bin/bash\n")
            file.write(f"#SBATCH --job-name={run_name}\n")
            file.write(f"#SBATCH --partition={self.partition}\n")
            file.write(f"#SBATCH --time={self.walltime_days}-{self.walltime_hours}:{self.walltime_mins}:{self.walltime_secs}\n")
            file.write(f"#SBATCH --nodes={self.nodes}\n")
            file.write(f"#SBATCH --ntasks={self.procs}\n")
            # file.write(f"#SBATCH --mail-user=erathnayake@sivananthanlabs.us\n")
            # file.write(f"#SBATCH --mail-type=ALL\n")
            file.write(f"mpirun -np {self.ntasks} pw.x -npool {self.npool} < {self.identifier}.in > {self.identifier}.out\n")
        os.rename(f"{self.identifier}.job", f"./{run_name}/{self.identifier}.job")

    def make_runs(self):
        """This method makes the runs. The inputs files are created in a method that handles the relevant file type
        These varibles have to be already set for this method to work:
        self.structure: Has to be set
            0 - The type of the structure is set from outside
            1 - An fcc cell that scales with the lattice constant
        """
        # Here the name version for the pseudopotentials is created since we have to have a form that can go on the folder names
        self.PP = ""
        self.specie = ""

        # Here we set the pseudo path according to what computer you would be running the code on
        self.set_pseudo_dir(self.job_handler)

        for key, val in self.pseudopotentials.items():
            self.PP = f"{val}-{self.PP}"
            self.specie = f"{key}-{self.specie}"
        if self.structure == None:
            print(f"make_runs: Warning! No structure set. Cannot create run.")
        else:
            for KE_cut_i in self.KE_cut:
                for a0_i in self.a0:
                    for k_i in self.k:
                        for R_i in self.R:  # This has been disabled for now
                            if R_i == None: 
                                R_name = f"{4*KE_cut_i}"
                            else:
                                R_name = R_i
                            run_name = f"{self.identifier}{self.d}Calc{self.equals}{self.calculator}{self.d}Struct{self.equals}{self.structure_type}{self.d}Specie{self.equals}{self.specie}{self.d}KE{self.equals}{KE_cut_i}{self.d}K{self.equals}{k_i}{self.d}R{self.equals}{R_name}{self.d}a{self.equals}{a0_i}{self.d}type{self.equals}{self.calculation}"
                            if self.structure == 0:
                                # cell has been set from outside
                                pass
                            elif self.structure == 1:
                                # An fcc cell that scales with the lattice constant
                                b = a0_i/2.0
                                self.atoms_object.set_cell([(0, b, b), (b, 0, b), (b, b, 0)], scale_atoms=True)
                            elif self.structure == 2:
                                # An simple cubic cell that scales with the lattice constant
                                b = a0_i
                                self.atoms_object.set_cell([(b, 0, 0), (0, b, 0), (0, 0, b)], scale_atoms=True)
                            else:
                                print("make_runs: Warning! Cannot set cell. Structrue not supported")
                            if os.path.exists(run_name):
                                shutil.rmtree(run_name)
                                print("make_runs: Warning! Path exists!! Overwriting")
                            os.mkdir(f"{run_name}")
                            self.write_QE_inputfile(run_name, KE_cut_i, R_i, a0_i, k_i)
                            self.all_runs_list.append(run_name)

                        # Creating jobs
                        # This if else block is pending deletion upon making a seperate method for slurm jobs.
                        if self.job_handler == "torque" or self.job_handler == "era_ubuntu" or self.job_handler == "era_pc" or self.job_handler == "sl_laptop":
                            pass
                            # we dont create job files for everyrun here for now since torque jobs will have specific script to run
                            # self.create_torque_job(run_name)
                        elif self.job_handler == "slurm":
                            self.create_slurm_job(run_name)
                        # elif self.job_handler == "era_pc":
                        #     self.create_era_pc_job(run_name)
                        else:
                            print(f"make_runs: Unrecognized job_handler! Job files not created")

        if self.job_handler == "torque":
            bat_file = open("rsyn_out.bat", "w+")
            bat_file.write(f'wsl rsync -avtuz -e "ssh -p 33301" ./ rathnayake@localhost:~/Run_files')
            bat_file.close()
            self.create_torque_job()
        elif self.job_handler == "era_ubuntu":
            bat_file = open("rsyn_out_eraubuntu.bat", "w+")
            bat_file.write(f'wsl rsync -avtuz -e ssh era@192.168.0.23 ./ era@192.168.0.23:~/Documents/Run_files')
            bat_file.close()
            self.create_job()
        elif self.job_handler == "era_pc" or self.job_handler == "sl_laptop":
            self.create_job()  # sine they are identical
        else:
            self.create_bash_file()

    def create_bash_file(self):
        """
        This script creates bash files so that you can run a batch of the runs that need to be done
        """
        print(f"create_bash_file: job_handler is set to: {self.job_handler}")
        bash_file = open(f"run_{self.identifier}.sh", "w+")
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
        elif self.job_handler == "slurm":
            bash_file.write(f'  sbatch *job\n')
        elif self.job_handler == "era_pc":
            bash_file.write(f'  . *job\n')
        bash_file.write(f'  cd ..\n')
        bash_file.write(f'done\n')
        bash_file.close()

class ReadOutfiles():
    """
    This method should read all out files of a given type (sesta/qe) and read the values like total energies
    identifier: list of strings
    """
    def __init__(self, *args, **kwargs):
        """
        All Kwargs should be set as strings
        MHP_base:
            This is the parameter that will be checked if the monkhorst pack grid is assymetric and the file folder has a list. 
            You will have to assign a char "x", "y", or "z" to this variable in order to let the program know what direction to chose to compare with others. 
            If you want to compare multiple you can just load another instance and load the requried direction.
            The value defaults to the "x" direction.
        """
        self.d = f"^"  # Here you can set the desired delimiter
        self.equals = ["+", "="]  # Here you can set the desired symbol for value assigner it can also be a list of all possible values

        # Gettings args here
            # No args to get right now

        # Getting kwargs here
        self.identifier       = kwargs.get("identifier", [])
        self.job_handler      = kwargs.get("job_handler", [])
        self.a0               = kwargs.get("a0", [])
        self.KE_cut           = kwargs.get("KE_cut", [])
        self.k                = kwargs.get("k", [])
        self.R                = kwargs.get("R", [])
        self.pseudopotentials = kwargs.get("pseudopotentials", [])
        self.pseudo_dir       = kwargs.get("pseudo_dir", [])
        self.calculator       = kwargs.get("calculator", [])
        self.structure_type   = kwargs.get("structure_type", [])
        self.xc               = kwargs.get("xc", [])
        self.calculation      = kwargs.get("calculation", [])
        self.species          = kwargs.get("species", [])
        self.high_verbosity   = kwargs.get("high_verbosity", False)
        self.MHP_base         = kwargs.get("MHP_base", "x")  # This is what you want to base you analysis on.

        # # Initializations
        self.atoms_objects = []  # Where all data read from file will be stored for a scf type file.
        self.atoms_bands_objects = []  #Where all data read from file will be stored for a bands type file.

    def read_folder_data(self, dir):
        """
        This method reads the out files from the requried directories
        """
        # print(f"list dir method {os.listdir()}")
        # for root, dirs, files in os.walk(os.getcwd(), topdown=False):
        #     for name in dirs:
        self.directory_list = os.listdir(dir)
        # print(f"Printing Directory list: {self.directory_list}")

        # Removing all the non relevant folders
        directoriestopop = []
        for dir in self.directory_list:
            if "^" not in dir:
                # print(f"read_folder_data: Ignoring folder or file (set to be removed): {dir}")
                directoriestopop.append(dir)
        # print(f"THis is the list of directories (before):\n{self.directory_list}")
        for x in directoriestopop:
            printable = self.directory_list.pop(self.directory_list.index(x))
            if self.high_verbosity == True:
                print(f"read_folder_data: Ignoring folder or file: {printable}")
        # print(f"THis is the list of directories(after):\n{self.directory_list}")

        self.folder_data = []
        # print("Right now we are in the read_folder_data method")  # For Debugging
        for dir in self.directory_list:
            # print(f"This is the dir: {dir}")  # For debugging
            try:
                x = dir.split("^")
                # print(x)
                run_parameters = {}
                run_parameters.update({"identifier":x[0]})  #Done seperately due to "identifier" not being present as a word in the folder name
                for i in range(1,len(x)):
                    y = x[i].split("=")
                    try:  # Any numerical values are tried to be converted float. If cannot just set as is.
                        run_parameters.update({y[0]:float(y[1])})
                    except:
                        run_parameters.update({y[0]:y[1]})

                # Adjusting for nonsymmetric monkhorst pack grids and basing the calling method to chose from "x", "y" or "z" directions in the Monkhorst Pack grid
                x = run_parameters["K"]
                if type(x) == str:  # This means the above except statement was executed
                    x = x.strip("]").strip("[").split(",")
                    if self.MHP_base == "x": x = int(x[0])
                    if self.MHP_base == "y": x = int(x[1])
                    if self.MHP_base == "z": x = int(x[2])
                run_parameters["K"] = x

                # Splitting up multiple values in Specie
                x = run_parameters["Specie"]
                x = x.split("-")
                x = [i for i in x if i != ""]
                run_parameters["Specie"] = x

                # # Splitting up multiple values in PP
                # x = run_parameters["PP"]
                # x = x.split("-")
                # x = [i for i in x if i != ""]
                # run_parameters["PP"] = x
                self.folder_data.append(run_parameters)

                # print(f"read_folder_data: Logging folder: {dir}")   #For debugging purposes
            except:
                print(f"read_folder_data: Warning!! Something wrong with {dir}. Cannot recognize patterns.")
            # print(run_parameters)  # For Debugging

    def make_required_folders_list(self):
        """
        This method reads the out file from the selected run files properties.
        Prerequisits:
            read_folder_data should be run before this
        """
        self.required_folders_list = []
        self.required_folder_data = []
        for folder in self.folder_data:
            if folder["identifier"] in self.identifier or self.identifier == []:
                print("inside Identifier")
                if folder["Calc"] in self.calculator or self.calculator == []:
                    # since you can mistakenly set a0 in strings lets try to convert them to floats
                    self.a0 = [float(x) for x in self.a0]
                    if folder["a"] in self.a0 or self.a0 == []:
                        if float(folder["R"]) in self.R or self.R ==[]:
                            if folder["Struct"] in self.structure_type or self.structure_type == []:
                                # Below commented lines for getting pseudos from teh file name but now not implemented
                                # count = 0
                                # for x in folder["PP"]:
                                #     if x in self.pseudopotentials or self.pseudopotentials == []:
                                #         count+=1
                                # if count == len(folder["PP"]):
                                if float(folder["KE"]) in self.KE_cut or self.KE_cut == []:
                                    if float(folder["K"]) in self.k or self.k == []:
                                        if folder["type"] in self.calculation or self.calculation == []:
                                            self.required_folders_list.append(self.directory_list[self.folder_data.index(folder)])
                                            self.required_folder_data.append(self.folder_data[self.folder_data.index(folder)])
        if self.high_verbosity == True:
            print("Loaded folders:")
            for folder in self.required_folders_list:
                print(folder)
    # def read_outfiles(self, directory, file_name):

    def read_outfiles(self, dir):
        """
        This method will handle the opening of all the files and then also put them together. 
        Inputs:
            dir (string): Can be "thesis", "sivalabs", or here"here"
            This will determine which folders to open. some folders default locations are hard coded.
        """

        # For debugging purposses below lines will be helpful
        # print("folders detected")
        # for x in self.directory_list:
        #     print(x)
        # print("folders to open")
        # for x in self.required_folders_list:
        #     print(x)

        from pathlib import Path

        cur_dir = Path(os.getcwd())
        runs_dir = cur_dir.parent.parent
        if dir == "thesis":
            mydir = Path(runs_dir, "Run_files_git", "Run_files")
        elif dir == "sivalabs":
            mydir = Path(runs_dir, "Run_files_SL", "Synced")
        elif dir== "here":
            mydir = Path(cur_dir)
        
        if self.high_verbosity:
            print(f"The Runs directory is: {mydir}")
        self.read_folder_data(mydir)
        self.make_required_folders_list()

        for x in range(0,len(self.required_folders_list)):
            path = os.path.join(mydir, self.required_folders_list[x], self.identifier[0])
            if self.high_verbosity:
                print(f"read_outfiles: Opening file: {path}.out")
            try:
                if self.folder_data[x]["Calc"].lower() == "qe":
                    file = ase.io.read(f"{path}.out", format = "espresso-out")
                    if self.calculation == "bands":
                        bands_file = ase.io.read(f"{path}.bands.out", format = "espresso-out")
                elif self.folder_data[x]["Calc"].lower() == "siesta":
                    file = ase.io.read(f"{path}.out", format = "espresso-out")
                self.atoms_objects.append(file)
                try:
                    self.atoms_bands_objects.append(bands_file)
                except:
                    if self.high_verbosity:
                        print(f"read_outfiles: Recognized as not a bands file. bands files not appeneded to atoms_bands_objects")
            except AssertionError as er:
                print(f"read_outfiles: ** Warning Fatal Error. 'espresso.py' in ASE is giving out an assertion error as below:")
                raise
            except:
                print(f"read_outfiles: ** Warning Fatal Error. Cannot read file. File might not be present or might not have finished Recommended to set parameters to specifically exclude this file.\n{path}.out\nIt might also be the bands file of the same name.\n{path}.bands.out")
        if self.atoms_bands_objects == []:
            # No bands files have been read
            self.data = list(zip(self.required_folders_list, self.required_folder_data, self.atoms_objects))
            if self.high_verbosity:
                print(f"read_outfiles: Sucessfully read scf files. Zipping done")
        else:
            # Trying to zip bands files here
            # print(self.atoms_bands_objects)
            self.data = list(zip(self.required_folders_list, self.required_folder_data, self.atoms_objects, self.atoms_bands_objects))
            if self.high_verbosity:
                print(f"read_outfiles: Sucessfully read bands files and scf files. Zipping done")

    def get_sorted_energies(self, sort_for = "KE"):
        """
        Returns two lists sorted according to sort_for
        """
        self.data = sorted(self.data, key=lambda x: x[1][sort_for])
        x = [self[1][sort_for] for self in self.data]
        E = [self[2].get_total_energy() for self in self.data]
        return x, E

    def get_band_path(self):
        if self.calculation == "scf":
            print("This is a scf calculation and therefore no band path.")
        elif self.calculation == "bands":
        # try:
            return self.atoms_bands_objects[0].cell.get_band_path()
        # except:
        #     print("Error returning band path")

def make_all_job_files(job_list = []):
    """
    This method makes all jobs run when executing a single file. 
    Warning: Not set to properly handle .bands files since .scf has to finish in order for the .bands files to run.
             Therefore functionality is only set for .scf.job files to be listed and according to how they finish you manually run the bands.job files
    """

    print("make_all_job_files: Printing all jobs onto a single file.")
    directory_list = os.listdir(os.getcwd())  # os.getcwd() might give different folders in different systems.
    with open("all_jobs.job", "w+") as file:
        file.write("#!/bin/bash\n\n")
        file.write("dos2unix *.job\n")
        # If there are no explicitly given jobs
        if len(job_list) == 0:
            for file2 in directory_list:
                # if "scf.job" in file2 or "bands.job" in file2:
                if "scf.job" in file2 or "relax.job" in file2:
                    file.write(f". {file2}\n")
        # if jobs are explicitly given
        else:
            for job in job_list:
                for file2 in directory_list:
                    # if job in file2 and (".scf.job" in file2 or ".bands.job" in file2):
                    if job in file2 and (".scf.job" in file2) or job in file2 and (".relax.job" in file2):
                        file.write(f". {file2}\n")
        file.write(f"\n")

if __name__ == "__main__":
    """
    This is used as an example as to how we use this file.
    """
    pass