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
        self.pseudopotentials = kwargs.get("pseudopotentials", {'Sn':'Sn_ONCV_PBE_FR-1.1.upf'})
        self.pseudo_dir       = kwargs.get("pseudo_dir", False)
        self.calculator       = kwargs.get("calculator", "QE")
        self.structure_type   = kwargs.get("structure_type", "bulk")
        self.xc               = kwargs.get("xc", "pbe")
        self.calculation      = kwargs.get("calculation", "scf")

        # self.R = kwargs.get("R", [300])

        # Quantum espresso inits
        self.ntasks          = kwargs.get("ntasks", 20)
        self.npool           = kwargs.get("npool", 1)
        self.espresso_inputs = {"pseudopotentials": self.pseudopotentials,
                                "calculation"     : self.calculation,
                                "lspinorb"        : kwargs.get("lspinorbit", False),
                                "noncolin"        : kwargs.get("noncolin", False),
                                # "ecutrho"         : KE_cut_i*4,
                                "occupations"     : kwargs.get("occupations",'smearing'),
                                "smearing"        : kwargs.get("smearing",'gaussian'),
                                "degauss"         : kwargs.get("degauss", 0.01),
                                "mixing_beta"     : kwargs.get("mixing_beta", 0.7),
                                "Title"           : kwargs.get("Title",'Sn'),
                                "prefix"          : kwargs.get("prefix",'Sn'),
                                "restart_mode"    : kwargs.get("restart_mode",'from_scratch'),
                                "disk_io"         : kwargs.get("disk_io",'default'),
                                "verbosity"       : kwargs.get("verbosity",'high'),
                                "lkpoint_dir"     : kwargs.get("lkpoint_dir", False),
                                "etot_conv_thr"   : kwargs.get("etot_conv_thr", 1.0e-6),
                                "forc_conv_thr"   : kwargs.get("forc_conv_thr", 1.0e-4),
                                "outdir"          : kwargs.get("outdir", './')
                                }

        if self.pseudo_dir != False:
            # Since if not set we want the value to be the default value
            # pseduo_dir is initialized in the base run init section
            self.espresso_inputs.update({"pseudo_dir" : self.pseudo_dir})

        # Here goes the job init stuff
        self.walltime_days   = kwargs.get("walltime_days", 2)
        self.walltime_mins   = kwargs.get("walltime_mins", 0)
        self.walltime_hours  = kwargs.get("walltime_hours", 2)
        self.walltime_secs   = kwargs.get("walltime_secs", 0)
        self.nodes           = kwargs.get("nodes", 2)
        self.procs           = kwargs.get("procs", 8)
        self.partition       = kwargs.get("partition", "cluster")

        # Other Initializations
            # For structure: An fcc cell that scales with the lattice constant = 1
        self.structure = kwargs.get("structure", 1)
        default = Atoms('Sn2', [(0, 0, 0), (0.25, 0.25, 0.25)],  pbc=True)
        self.atoms_object = kwargs.get("atoms_object", default)
        self.all_runs_list = []

    def set_pseudo_dir(self, machine):
        """
        Sets the pseudo_dir according to the machine
        """
        pseudo_database_path = {"cluster":"/usr/local/share/espresso/pseudo",
                        "carbon":"../PseudopotentialDatabase",
                        "siva_labs_wsl":"/mnt/c/Users/Eranjan/Desktop/PseudopotentialDatabase",
                        "home_wsl":"/mnt/c/Users/Eranjan/Desktop/PseudopotentialDatabase"
                        }
        self.espresso_inputs.update({"pseudo_dir" : pseudo_database_path[machine]})

    def set_pseudopotentials(self, pseudos):
        """
        Sets the pseudopotentials
        """
        self.pseudopotentials = pseudos
        self.espresso_inputs.update({"pseudopotentials": pseudos})

    def write_QE_inputfile(self, run_name, KE_cut_i, a0_i, k_i):
        """
        This method creates Quantum espresso input files
        """
        # Here we update aditional run specific stuff
        self.espresso_inputs.update({"label" : f"{run_name}"})
        self.espresso_inputs.update({"kpts" : (k_i, k_i, k_i)})
        self.espresso_inputs.update({"ecutwfc" : KE_cut_i})
        ase.io.write(f"{self.identifier}.in", self.atoms_object, format = "espresso-in", **self.espresso_inputs)

    def write_SIESTA_inputfile(self, run_name, KE_cut_i, a0_i, k_i):
        """
        This method creates SIESTA input files
        """
        pass

    def get_number_of_calculations(self):
        return (len(self.KE_cut)*len(self.a0)*len(self.k))

    def create_torque_job(self):
        with open (f"{self.identifier}.job", "w+") as file_torque:
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
            file_torque.write(f"    mpirun pw.x < {self.identifier}.in > {self.identifier}.out\n")
            file_torque.write(f"END_JOB_SCRIPT\n")
            file_torque.write(f"\n")
            file_torque.write(f"done <<'END_TASKLIST'\n")
            file_torque.write(f"    # Single quoting the limit string 'EOT' will pass strings without shell variable and execution expansion.\n")
            file_torque.write(f"    # Comments and empty line are fine because we explicitly skip them.\n\n")
            for run in self.all_runs_list:
                file_torque.write(f"    {run}\n")
            file_torque.write(f"END_TASKLIST\n")
        # os.rename(f"{self.identifier}.job", f"./{run_name}/{self.identifier}.job")

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
                        run_name = f"{self.identifier}{self.d}Calc{self.equals}{self.calculator}{self.d}Struct{self.equals}{self.structure_type}{self.d}Specie{self.equals}{self.specie}{self.d}KE{self.equals}{KE_cut_i}{self.d}K{self.equals}{k_i}{self.d}R{self.equals}{R_i}{self.d}a{self.equals}{a0_i}{self.d}PP{self.equals}{self.PP}{self.d}type{self.equals}{self.calculation}"
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
                            # we do nothing here for now since torque jobs will have specific script to run
                            pass
                            # self.create_torque_job(run_name)
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

class ReadOutfiles():
    """This method should read all out files of a given type (sesta/qe) and read the vlaues like total energies"""
    def __init__(self, *args, **kwargs):
        """All Kwargs should be set as strings"""
        self.d = f"^"  # Here you can set the desired delimiter
        self.equals = ["+", "="]  # Here you can set the desired symbol for value assigner it can also be a list of all possible values

        # Gettings args here

        # Getting kwargs here
        self.identifier       = kwargs.get("identifier", [])
        self.job_handler      = kwargs.get("job_handler", [])
        self.a0               = kwargs.get("a0", [])
        self.KE_cut           = kwargs.get("KE_cut", [])
        self.k                = kwargs.get("k", [])
        self.pseudopotentials = kwargs.get("pseudopotentials", [])
        self.pseudo_dir       = kwargs.get("pseudo_dir", [])
        self.calculator       = kwargs.get("calculator", [])
        self.structure_type   = kwargs.get("structure_type", [])
        self.xc               = kwargs.get("xc", [])
        self.calculation      = kwargs.get("calculation", [])
        self.species          = kwargs.get("species", [])

        # # Initializations
        self.atoms_objects = []

    def read_folder_data(self, dir):
        """
        This method reads the out files from the requried directories
        """
        # print(f"list dir method {os.listdir()}")
        # for root, dirs, files in os.walk(os.getcwd(), topdown=False):
        #     for name in dirs:
        self.directory_list = os.listdir(dir)
        self.folder_data = []
        # print("Right now we are in the read_folder_data method")  # For Debugging
        for dir in self.directory_list:
            print(f"This is the dir: {dir}")  # For debugging
            try:
                x = dir.split("^")
                run_parameters = {}
                run_parameters.update({"identifier":x[0]})  #Done seperately due to "identifier" not being present as a word in the folder name
                for i in range(1,len(x)):
                    y = x[i].split("=")
                    run_parameters.update({y[0]:y[1]})
                # print(run_parameters)
                # Splitting up multiple values in Specie
                x = run_parameters["Specie"]
                x = x.split("-")
                x = [i for i in x if i != ""]
                run_parameters["Specie"] = x

                # Splitting up multiple values in PP
                x = run_parameters["PP"]
                x = x.split("-")
                x = [i for i in x if i != ""]
                run_parameters["PP"] = x
                self.folder_data.append(run_parameters)
                # print(f"read_folder_data: Logging folder: {dir}")   #For debugging purposes
            except:
                print(f"read_folder_data: Ignoring folder/file: {self.directory_list.pop(self.directory_list.index(dir))}")
            # self.directory_list.pop(self.directory_list.index("make_ligands.py"))
        # print(self.directory_list)
        # print(self.folder_data)

            ## here we need to flush out other directries that are not run directories

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
                # print(self.identifier)  # for debugging purposes
                if folder["Calc"] in self.calculator or self.calculator == []:
                    if folder["a"] in self.a0 or self.a0 == []:
                        if folder["Struct"] in self.structure_type or self.structure_type == []:
                            count = 0
                            for x in folder["PP"]:
                                if x in self.pseudopotentials or self.pseudopotentials == []:
                                    count+=1
                            # print(count)
                            if count == len(folder["PP"]):
                                if folder["KE"] in self.KE_cut or self.KE_cut == []:
                                    if folder["K"] in self.k or self.k == []:
                                        if folder["type"] in self.calculation or self.calculation == []:
                                            self.required_folders_list.append(self.directory_list[self.folder_data.index(folder)])
                                            self.required_folder_data.append(self.folder_data[self.folder_data.index(folder)])

    # def read_outfiles(self, directory, file_name):

    def read_outfiles(self, dir):
        """
        This method reads the out file from a single run / single folder 
        """

        # For debugging purposses below lines will be helpful
        # print("folders detected")
        # for x in self.directory_list:
        #     print(x)
        # print("folders to open")
        # for x in self.required_folders_list:
        #     print(x)

        from pathlib import Path
        # from ebk.convergence import E_cut_Optimize

        cur_dir = Path(os.getcwd())
        runs_dir = cur_dir.parent.parent
        if dir == "thesis":
            mydir = Path(runs_dir, "Run_files")
        elif dir == "sivalabs":
            mydir = Path(runs_dir, "Run_files_SL", "Synced")
        print(f"The Runs directory is: {mydir}")

        self.read_folder_data(mydir)
        self.make_required_folders_list()

        for x in range(0,len(self.required_folders_list)):
            path = os.path.join(mydir, self.required_folders_list[x], self.identifier[0])
            print(f"Opening file: {path}.out")
            file = ase.io.read(f"{path}.out", format = "espresso-out")
            self.atoms_objects.append(file)
        print(self.atoms_objects)

            # try:
            #     file = ase.io.read(f"{path}.out", format = "espresso-out")
            #     self.atoms_objects.append(file)
            # except:
            #     print(f"read_folder_data: Ignoring folder/file: { self.required_folders_list.pop(x)}")
            #     print(f"Probably no out file")

if __name__ == "__main__":
    """This is used as an example as to how we use this file."""
    pass
