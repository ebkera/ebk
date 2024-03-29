"""
This file creates:
SIESTAinput files
QE input files
bash scripts for running on torque.
bash scripts for running on local machines
This can be used to run multiple jobs for example
As of now this code can only do runs with similar job scheduler parameters
"""

import os, sys
from ase.build import bulk
from ase import Atoms
import ase.io
from ase.io import read, write
from ebk.QE import QErunfilecreator  # So that we can see how we did it last time
import shutil
from ebk.SIESTA.SIESTAcontrol import Generatefdf

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
            "k" (list of lists with length 3 or list of ints): The k grid
            "k_nscf" (list of length 3 or single int): The k grid that you want to preserve for nscf calculations by default it will be the same as the relevant scf calculation
            "pseudopotentials" (string):
            "atoms_object" (atoms object): This should be without setting the cell since that will be done with every iteration
            "structure" (int): vlaue will determine the cell
                1: fcc structure with a0 as the lattice constant
        """
        # Here we have some variables for the script
        self.d = f"^"  # Here you can set the desired delimiter
        self.dm = f"+"  # Minor delimiter. This is the minor delimiter where you can delimit other things (used in type (self.calculation) eg: ^type=scf+dos+nscf+pdos)
        self.equals = f"="  # Here you can set the desired symbol for value assigner\
        self.files_to_move = []
        self.files_to_copy = []
        self.files_to_move_to_runfolder = []
        self.files_to_copy_to_runfolder = []

        # Gettings args here
        # No args to get currently

        # Setting kwargs here
        # Base Run inits
        self.identifier       = kwargs.get("identifier", "run")
        self.SystemLabel       = kwargs.get("SystemLabel", "run")
        if "identifier" not in kwargs.keys(): self.identifier       = kwargs.get("SystemLabel", "run")  # This is if a siesta input is given
        self.job_handler      = kwargs.get("job_handler", "")
        self.a0               = kwargs.get("a0", [6.47])
        self.KE_cut           = kwargs.get("KE_cut", [20])
        if "PAO_EnergyShift" in kwargs.keys(): self.KE_cut       = [kwargs.get("PAO_EnergyShift", 0.001)]  # This is if a siesta input is given
        self.k                = kwargs.get("k", [2])
        self.k_nscf           = kwargs.get("k_nscf", self.k[0])
        self.pseudopotentials = kwargs.get("pseudopotentials", {'Sn':'Sn_ONCV_PBE_FR-1.1.upf'})
        self.pseudo_dir       = kwargs.get("pseudo_dir", False)
        self.calculator       = kwargs.get("calculator", "QE").upper()
        self.structure_type   = kwargs.get("structure_type", "bulk")
        if "fdf_type" in kwargs.keys(): self.structure_type       = kwargs.get("fdf_type", "bulk")  # This is if a siesta input is given
        self.xc               = kwargs.get("xc", "pbe")
        self.calculation      = kwargs.get("calculation", "scf")
        self.path             = kwargs.get("path", "GXWLGKL")
        self.density          = kwargs.get("density", 20)
        self.k_path           = {"path":self.path, "density": self.density}
        self.R                = kwargs.get("R", [None])
        self.base_folder      = kwargs.get("base_folder", "Runs")
        if "occupations" in kwargs: self.occupations      = kwargs.get("occupations",'smearing')
        if "occupations_nscf" in kwargs: self.occupations_nscf = kwargs.get("occupations_nscf",'tetrahedra')
        # For DOS PDOS and LDOS calculations. This does not use ASE.
        self.DOS_EMIN       = kwargs.get("DOS_EMIN",-20)
        self.DOS_EMAX       = kwargs.get("DOS_EMAX",20)
        self.DOS_DeltaE     = kwargs.get("DOS_DeltaE",0.1)
        self.DOS_degauss    = kwargs.get("DOS_degauss",0.01)
        self.PDOS_EMIN       = kwargs.get("PDOS_EMIN", self.DOS_EMIN)
        self.PDOS_EMAX       = kwargs.get("PDOS_EMAX", self.DOS_EMAX)
        self.PDOS_DeltaE     = kwargs.get("PDOS_DeltaE", self.DOS_DeltaE)
        self.PDOS_degauss    = kwargs.get("PDOS_degauss", self.DOS_degauss)
        self.PDOS_required_projections = kwargs.get("PDOS_required_projections", [["Hg","all"],["Hg","s"],["Hg","p"],["Hg","d"],["Te","all"],["Te","s"],["Te","p"],["Te","d"],["H","s"],["S","all"],["S","s"],["S","p"],["C","all"],["C","s"],["C","p"]])
        self.LDOS_EMIN       = kwargs.get("PDOS_EMIN", self.DOS_EMIN)
        self.LDOS_EMAX       = kwargs.get("PDOS_EMAX", self.DOS_EMAX)
        self.LDOS_DeltaE     = kwargs.get("PDOS_DeltaE", self.DOS_DeltaE)
        self.LDOS_degauss    = kwargs.get("PDOS_degauss", self.DOS_degauss)
        self.fft_grid        = kwargs.get("fft_grid", [40, 40, 432])
        self.n_proj_boxes    = kwargs.get("n_proj_boxes", self.fft_grid[2])
        self.job_name        = kwargs.get("job_name", self.identifier)
        self.email_recipients= kwargs.get("email_recipients", ["ebk_era@hotmail.com"])  # List of strings
        self.send_mail       = kwargs.get("send_mail", True)

        # Quantum espresso inits some other inits that need to be only set if explicitly given can be found below this.
        self.espresso_inputs = {"pseudopotentials": self.pseudopotentials,
                                "calculation"     : self.calculation,
                                "lspinorb"        : kwargs.get("lspinorb", False),
                                "noncolin"        : kwargs.get("noncolin", False),
                                "diagonalization" : kwargs.get("diagonalization",'david'),
                                "mixing_beta"     : kwargs.get("mixing_beta", 0.2),
                                "Title"           : kwargs.get("Title",'Sn'),
                                "prefix"          : kwargs.get("prefix",'Sn'),
                                "verbosity"       : kwargs.get("verbosity",'high'),
                                "wf_collect"      : kwargs.get("wf_collect", False),
                                # "disk_io"         : kwargs.get("disk_io",'low'),
                                # "disk_io_nscf"    : kwargs.get("disk_io_nscf",'low'),
                                "etot_conv_thr"   : kwargs.get("etot_conv_thr", 1.0e-5),
                                "forc_conv_thr"   : kwargs.get("forc_conv_thr", 1.0e-5),
                                "outdir"          : kwargs.get("outdir", './'),
                                "path"            : self.path,
                                "density"         : self.density,  # This is an ASE command for input files for Quantum Espresso
                                "electron_maxstep": kwargs.get("electron_maxstep", 1000),
                                "nstep"           : kwargs.get("nstep", 1000),
                                "mixing_mode"     : kwargs.get("mixing_mode", "plain"),
                                "cell_factor"     : kwargs.get("cell_factor", 1.2),
                                "cell_dofree"     : kwargs.get("cell_dofree", 'z'),
                                }

        if self.calculator == "SIESTA":
            self.SIESTA_inputs = kwargs 

        # Here are all initializations of the self.espresso_inputs variable that should be set only if explicitly given by user
        # Here as you can see the default values sef for the kwargs will never get used. They are there as a guide to what you can use when you might need to use them
        # Also another this is the values here does not necessarily equal the defalt values for the variables in Quantum Espresso
        if "restart_mode" in kwargs:
            self.espresso_inputs.update({"restart_mode"    : kwargs.get("restart_mode",'from_scratch')})
        if "lkpoint_dir" in kwargs:
            self.espresso_inputs.update({"lkpoint_dir"     : kwargs.get("lkpoint_dir", False)})
        if "nbnd" in kwargs:
            self.espresso_inputs.update({"nbnd"            : kwargs.get("nbnd", 40)})
        if "degauss" in kwargs:
            self.espresso_inputs.update({"degauss"         : kwargs.get("degauss", 0.01)})
        if "occupations" in kwargs:
            self.espresso_inputs.update({"occupations"     : kwargs.get("occupations", "smearing")})
        if "smearing" in kwargs:
            self.espresso_inputs.update({"smearing"        : kwargs.get("smearing", "gaussian")})
        if "cell_dynamics" in kwargs:
            self.espresso_inputs.update({"cell_dynamics"   : kwargs.get("cell_dynamics", "none")})         

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

        # General Prallelization
        self.nodes           = kwargs.get("nodes", 2)  # the number of cores
        self.procs           = kwargs.get("procs", 8)  # number of processesors per core
        self.partition       = kwargs.get("partition", "cluster")  # The partition that the job will run on
        self.ntasks          = kwargs.get("ntasks", 20)  # number of threads in total This is a slurm command

        # Quantum espresso parallelization
        self.nimage          = kwargs.get("nimage", 1)
        self.npools          = kwargs.get("npools", 1)  # The number of pools of k points per proc (has to be an integer). This is a Quantum espresso parameter and will only work with QE. 
        self.nband           = kwargs.get("nband", 1)
        self.ntg             = kwargs.get("ntg", 1)

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
                        "sl_laptop": "/mnt/c/Users/erathnayake/Desktop/PseudopotentialDatabase",
                        # "slurm": "/home/erathnayake/Synced/PseudopotentialDatabase",
                        }
        if machine in pseudo_database_path:
            # This means that if machine is not in above dict the pseudo path will not be updated.
            # This means that the system specific default folder will be chosen
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
            # print(f"inside lists")
            # print(k_i)
            self.espresso_inputs.update({"kpts" : (k_i[0], k_i[1], k_i[2])})
        else:
            self.espresso_inputs.update({"kpts" : (k_i, k_i, k_i)})
        if R_i != None: self.espresso_inputs.update({"ecutrho" : R_i})

        if "scf" in self.calculation or "bands" in self.calculation or "nscf" in self.calculation:
            # First we deal with the scf run.
            # The relax runs will also be saved with the .scf.out extension
            self.espresso_inputs.update({"calculation" : "scf"})
            ase.io.write(f"{self.identifier}.scf.in", self.atoms_object, format = "espresso-in", **self.espresso_inputs)
            os.rename(f"{self.identifier}.scf.in", f"./{self.base_folder}/{run_name}/{self.identifier}.scf.in")
        if "relax" in self.calculation and "vc" not in self.calculation:
            # This is a relax calculation.
            # The relax runs will also be saved with the .scf.out extension
            self.espresso_inputs.update({"calculation" : "relax"})
            ase.io.write(f"{self.identifier}.scf.in", self.atoms_object, format = "espresso-in", **self.espresso_inputs)
            os.rename(f"{self.identifier}.scf.in", f"./{self.base_folder}/{run_name}/{self.identifier}.scf.in")
        if "relax" in self.calculation and "vc" in self.calculation:
            # This is a vc relax calculation.
            self.espresso_inputs.update({"calculation" : "vc-relax"})
            ase.io.write(f"{self.identifier}.scf.in", self.atoms_object, format = "espresso-in", **self.espresso_inputs)
            os.rename(f"{self.identifier}.scf.in", f"./{self.base_folder}/{run_name}/{self.identifier}.scf.in")
        if "bands" in self.calculation:
            # Then with the bands file
            self.espresso_inputs.update({"calculation" : "bands"})
            self.espresso_inputs.update({"kpts" : self.k_path})
            ase.io.write(f"{self.identifier}.bands.in", self.atoms_object, format = "espresso-in", **self.espresso_inputs)
            os.rename(f"{self.identifier}.bands.in", f"./{self.base_folder}/{run_name}/{self.identifier}.bands.in")
        if "nscf" in self.calculation:
            # First we deal with the scf run. This will overwrite any scf runs produced during bands file creation.
            # if os.path.exists(f"./{self.base_folder}/{run_name}/{self.identifier}.scf.in"):
            #     os.remove(f"./{self.base_folder}/{run_name}/{self.identifier}.scf.in")
            # self.espresso_inputs.update({"calculation" : "scf"})

            # Then here we take care of the nscf run
            self.espresso_inputs.update({"calculation" : "nscf"})
            if type(self.k_nscf) == list:
                self.espresso_inputs.update({"kpts" : (self.k_nscf[0], self.k_nscf[1], self.k_nscf[2])})
            else:
                self.espresso_inputs.update({"kpts" : (self.k_nscf, self.k_nscf, self.k_nscf)})
            if self.occupations_nscf == "tetrahedra":
                try:
                    # A try here since smearing might not be present.
                    del(self.espresso_inputs["smearing"])
                except: pass
                try:
                    # A try here since smearing might not be present.
                    del(self.espresso_inputs["degauss"])
                except: pass
            self.espresso_inputs.update({"occupations" : "tetrahedra"})
            if "disk_io_nscf" in self.espresso_inputs.keys():
                self.espresso_inputs["disk_io"] = self.espresso_inputs["disk_io_nscf"]
            ase.io.write(f"{self.identifier}.nscf.in", self.atoms_object, format = "espresso-in", **self.espresso_inputs)
            # Putting files into folders
            os.rename(f"{self.identifier}.nscf.in", f"./{self.base_folder}/{run_name}/{self.identifier}.nscf.in")
        if f"{self.dm}dos":
            # Here we have to be careful since if we just use "dos" in the if statement we will get all teh pdos and ldos stuff as well
            # We are not using ASE to write this file. It is simple and that is one reason
            with open(f"{self.identifier}.dos.in", "w+") as file:
                file.write(f"&DOS\n")
                file.write(f"  prefix  = '{self.espresso_inputs['prefix']}'\n")
                file.write(f"  fildos  = '{self.espresso_inputs['prefix']}.dos.dat'\n")
                file.write(f"  Emin    = {self.DOS_EMIN}\n")
                file.write(f"  Emax    = {self.DOS_EMAX}\n")
                file.write(f"  DeltaE  = {self.DOS_DeltaE}\n")
                file.write(f"  degauss = {self.DOS_degauss}\n")
                file.write(f"/\n")
            os.rename(f"{self.identifier}.dos.in", f"./{self.base_folder}/{run_name}/{self.identifier}.dos.in")  
        if "pdos" in self.calculation:
            # We are not using ASE to write this file. It is simple and that is one reason
            with open(f"{self.identifier}.pdos.in", "w+") as file:
                file.write(f"&projwfc\n")
                file.write(f"  prefix  = '{self.espresso_inputs['prefix']}'\n")
                file.write(f"  filpdos = '{self.espresso_inputs['prefix']}.pdos.dat'\n")
                file.write(f"  Emin    = {self.PDOS_EMIN}\n")
                file.write(f"  Emax    = {self.PDOS_EMAX}\n")
                file.write(f"  DeltaE  = {self.PDOS_DeltaE}\n")
                file.write(f"  degauss = {self.PDOS_degauss}\n")
                file.write(f"/\n")
            os.rename(f"{self.identifier}.pdos.in", f"./{self.base_folder}/{run_name}/{self.identifier}.pdos.in")
        if "ldos" in self.calculation:
            # We are not using ASE to write this file. It is simple and that is one reason
            with open(f"{self.identifier}.ldos.in", "w+") as file:
                file.write(f"&projwfc\n")
                file.write(f"  prefix       = '{self.espresso_inputs['prefix']}',\n")
                # file.write(f"  filpdos = '{self.espresso_inputs['prefix']}'\n")
                # file.write(f"  Emin    = {self.PDOS_EMIN}\n")
                # file.write(f"  Emax    = {self.PDOS_EMAX}\n")
                file.write(f"  DeltaE       = {self.LDOS_DeltaE},\n")
                file.write(f"  degauss      = {self.LDOS_degauss},\n")
                file.write(f"  tdosinboxes  = .true.,\n")
                file.write(f"  plotboxes    = .true.,\n")
                file.write(f"  n_proj_boxes = {self.n_proj_boxes},\n")
                for i in range(1, self.n_proj_boxes+1):
                    file.write(f"  irmin(1,{i})=1, irmax(1,{i})={self.fft_grid[0]}, irmin(2,{i})=1, irmax(2,{i})={self.fft_grid[1]}, irmin(3,{i})={i}, irmax(3,{i})={i},\n")
                file.write(f"/\n")
            os.rename(f"{self.identifier}.ldos.in", f"./{self.base_folder}/{run_name}/{self.identifier}.ldos.in")  

    def write_SIESTA_inputfile(self, run_name):
        """
        This method creates SIESTA input files
        """
        fdf = Generatefdf(**self.SIESTA_inputs)
        fdf.write()
        os.rename(f"{self.SIESTA_inputs['SystemLabel']}.fdf", f"./{self.base_folder}/{run_name}/{self.SIESTA_inputs['SystemLabel']}.fdf")
        if "Write.Denchar" in self.SIESTA_inputs.keys():
            if self.SIESTA_inputs["Write.Denchar"]:
                fdf.write_denchar()
                os.rename(f"{self.SIESTA_inputs['SystemLabel']}.Denchar.fdf", f"./{self.base_folder}/{run_name}/{self.SIESTA_inputs['SystemLabel']}.Denchar.fdf")
        
        # Copying PP files for siesta
        from ebk import get_machine_paths
        paths = get_machine_paths()
        # os.rename(f"{self.identifier}.scf.in", f"./{self.base_folder}/{run_name}/{self.identifier}.scf.in")
        for x in self.SIESTA_inputs["Species"]:
            y = x.split("_")[0]
            # Here we copy the required PP files. Ideally this can be pulled from the net but we need to work offline. Use git to control the content and use ebk.get_achine_paths()
            if y == "Sn" or y == "Pd" or y=="Se":
                if self.SIESTA_inputs['XC_Functional'] == "GGA" or self.SIESTA_inputs['XC_Functional'] == "VDW": xx = "GGA"
                else: xx = "LDA"
                shutil.copy(f"{paths['pps']}/_{xx}_rivero/{y}.psf", f'./{self.base_folder}/{run_name}/{x}.psf')
            else:
                shutil.copy(f"{paths['pps']}/{self.SIESTA_inputs['XC_Functional']}_PSF/{y}.psf", f'./{self.base_folder}/{run_name}/{x}.psf')

    def move_files_to_run_folder(self, run_name):
        """
        This method addf functionality to add any other file into a runfolder. This is will enable adding files generated seperately into the folder.
        Examples are pp files or maybe xyz or fdf files.
        Functionality:
            Copies files that are already in the main folder where this script is run into the base_folder/{run_name}/
        """
        if len(self.files_to_copy) == 0 and len(self.files_to_move) == 0: return
        else:
            for file_name in self.files_to_move:
                try:
                    shutil.move(file_name, f'./{self.base_folder}/{run_name}/{file_name}')
                except:
                    print(f"could not move file (might not be present/created): {file_name}")
            for file_name in self.files_to_copy:
                try:
                    shutil.copy(file_name, f'./{self.base_folder}/{run_name}/{file_name}')
                except:
                    print(f"could not coply file (might not be present/created): {file_name}")

        if len(self.files_to_copy_to_runfolder) == 0 and len(self.files_to_move_to_runfolder) == 0: return
        else:
            for file_name in self.files_to_move_to_runfolder:
                try:
                    shutil.move(file_name, f'./{self.base_folder}/{file_name}')
                except:
                    print(f"could not move file (might not be present/created): {file_name}")
            for file_name in self.files_to_copy_to_runfolder:
                try:
                    shutil.copy(file_name, f'./{self.base_folder}/{file_name}')
                except:
                    print(f"could not coply file (might not be present/created): {file_name}")

    def get_number_of_calculations(self):
        return (len(self.KE_cut)*len(self.a0)*len(self.k)*len(self.R))

    def create_job(self):
        """
        This is the most generic job creator. Will create all jobs
        """
        def create_runline(self):
            """On going project for bringing in all the line creators together. Not implemented"""
            if "bands" in self.calculation:
                file_torque.write(f"    mpirun pw.x < {self.identifier}.bands.in > {self.identifier}.bands.out\n")
                file_torque.write(f"    now=$(date)\n")
                file_torque.write(f'    echo "$now: $dir : completed bands" >> ../all_jobs.log\n')
            if "nscf" in self.calculation:
                file_torque.write(f"    mpirun -np {self.ntasks} pw.x < {self.identifier}.nscf.in > {self.identifier}.nscf.out\n")
                file_torque.write(f"    now=$(date)\n")
                file_torque.write(f'    echo "$now: $dir : completed nscf" >> ../all_jobs.log\n')

        with open (f"{self.base_folder}/{self.identifier}.{self.calculation}.job", "w+") as file_torque:
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
            # file_torque.write("    job_name=${dir#*KE=}\n")
            # file_torque.write("    job_name=${job_name%^a*}\n")
            file_torque.write(f"    job_name={self.job_name}")
            file_torque.write(f"    # Use another here-doc to read the job script, to obviate individual files\n")
            file_torque.write(f"    # in each data directory.\n")
            file_torque.write(f"\n")
            file_torque.write(f"    # A here-doc with leading '-' will get leading TABs removed.  (unwise to\n")
            file_torque.write(f"    # use, however, if your editor munges TABs.)\n")
            file_torque.write(f"\n")
            file_torque.write(f"    # Double-quote $dir to avoid parameter substitution.  (This is not strictly\n")
            file_torque.write(f"    # necessary here, though, because the chars '^+' are not special -- in the\n")
            file_torque.write(f"    # circumstances used here.)\n")
            file_torque.write(f"    #!/bin/bash\n")
            # file_torque.write(f"    now=$(date)\n")
            # file_torque.write(f'    echo "$now : $dir : Starting" >> all_jobs.log\n')
            if self.job_handler == "torque":
                file_torque.write(f'    qsub -w "$PWD/$dir" -N "$job_name" <<-END_JOB_SCRIPT\n')
                file_torque.write(f"    #  Basics: Number of nodes, processors per node (ppn), and walltime (hhh:mm:ss)\n")
                file_torque.write(f"    #PBS -l nodes={self.nodes}:ppn={self.procs}\n")
                file_torque.write(f"    #PBS -l walltime={self.walltime_hours}:{self.walltime_mins}:{self.walltime_secs}\n")
                # file.write(f"    #PBS -N {run_name}\n")
                file_torque.write(f"    #PBS -A cnm77824\n")
                file_torque.write(f"    #  File names for stdout and stderr.  If not set here, the defaults\n")
                file_torque.write(f"    #  are <JOBNAME>.o<JOBNUM> and <JOBNAME>.e<JOBNUM>\n")
                # file.write(f"    #PBS -o job.out\n")
                file_torque.write(f"    #PBS -e $PWD/$dir/job.err\n")
                file_torque.write(f"    # Join the output and error streams into the standard error file\n")
                file_torque.write(f"    #PBS -j eo\n")
                file_torque.write(f"    # Send mail at begin, end, abort, or never (b, e, a, n). Default is 'a'.\n")
                file_torque.write(f"    #PBS -m bea\n")
                file_torque.write("\n")
                # file_torque.write(f"    # change into the directory where qsub will be executed\n")
                file_torque.write(f"    cd \$PBS_O_WORKDIR\n\n")
                file_torque.write(f"    # use a per-job lineup of modules; stern\n")
                file_torque.write(f"    module purge\n")
                if self.calculator == "QE":
                    file_torque.write(f"    module load intel\n")
                    file_torque.write(f"    module load openmpi/1.10/intel-17\n")
                    file_torque.write(f"    module load quantum-espresso/5.4/openmpi-1.10\n")
                if self.calculator == "SIESTA":
                    file_torque.write(f"    module load intel/18\n")
                    # file_torque.write(f"    module load impi/5/\n")
                    file_torque.write(f"    module load impi/\n")  # Loads the impi18 module
                    # file_torque.write(f"    module load fftw3/3.3/impi-5\n")  # this is what the 4.1-b2-2 uses
                    file_torque.write(f"    module load fftw3/3.3/intel/3.3.2-1\n")  # this is what the 4.1-b2-2 uses
                    # file_torque.write(f"    module load siesta\n")
                file_torque.write(f"    module list\n\n")
                file_torque.write(f'    echo "PBS_O_WORKDIR: $PBS_O_WORKDIR"\n')
                file_torque.write(f'    echo "************** Starting Calculation ***************"\n')
                file_torque.write(f'    echo "Dir: $dir"\n')
                # file.write(f"    # start MPI job over default interconnect; count allocated cores on the fly.\n")
                # file.write(f"    mpirun -machinefile  $PBS_NODEFILE -np $PBS_NP pw.x < {run_name}.in > {run_name}.out\n")
                # file_torque.write(f"    now=$(date)\n")
                #Here Npools vs npool
                file_torque.write(f'    date\n')
                if self.calculator == "QE":
                    file_torque.write(f'    echo "Starting scf"\n')
                    file_torque.write(f"    mpirun pw.x -npools {self.npools} -ntg {self.ntg} -in {self.identifier}.scf.in > {self.identifier}.scf.out\n")
                    file_torque.write(f'    date\n')
                    file_torque.write(f'    echo "Completed scf"\n')
                    file_torque.write(f'    rm *wfc*\n')
                    file_torque.write(f'    echo "Removed wavefunction files"\n')
                    if "bands" in self.calculation:
                        file_torque.write(f'    date\n')
                        file_torque.write(f'    echo "Starting bands"\n')
                        file_torque.write(f"    mpirun pw.x -npools {self.npools} -ntg {self.ntg} -in {self.identifier}.bands.in > {self.identifier}.bands.out\n")
                        file_torque.write(f'    date\n')
                        file_torque.write(f'    echo "Completed bands"\n')
                    if "nscf" in self.calculation:
                        file_torque.write(f'    date\n')
                        file_torque.write(f'    echo "Starting nscf"\n')
                        file_torque.write(f"    mpirun pw.x -npools {self.npools} -ntg {self.ntg} -in {self.identifier}.nscf.in > {self.identifier}.nscf.out\n")
                        file_torque.write(f'    date\n')
                        file_torque.write(f'    echo "Completed nscf"\n')
                    if "pdos" in self.calculation:
                        file_torque.write(f'    echo "Calculationg PDOS"\n')
                        file_torque.write(f"    mpirun projwfc.x < {self.identifier}.pdos.in > {self.identifier}.pdos.out\n")
                        file_torque.write(f'    echo "Calculating PDOS components"\n')
                        for x in self.PDOS_required_projections:
                            if x[1] == "all":
                                file_torque.write(f"    sumpdos.x *\({x[0]}\)* > {self.identifier}.{x[0]}_all.PDOS\n")
                            else:
                                file_torque.write(f"    sumpdos.x *\({x[0]}\)*\({x[1]}\) > {self.identifier}.{x[0]}_{x[1]}.PDOS\n")
                    if f"{self.dm}dos" in self.calculation:
                        file_torque.write(f'    echo "Calculating DOS"\n')
                        file_torque.write(f"    mpirun dos.x < {self.identifier}.dos.in > {self.identifier}.dos.out\n")
                    if "ldos" in self.calculation:
                        file_torque.write(f'    echo "Calculating LDOS"\n')
                        file_torque.write(f"    mpirun projwfc.x < {self.identifier}.ldos.in > {self.identifier}.ldos.out\n")
            # file.write(f'rm *wfc*\n')
                if self.calculator == "SIESTA":
                    # file_torque.write(f"    ln -s ~/SIESTA_compile/siesta-master/Obj/siesta siesta_v0\n")
                    file_torque.write(f"    mpirun ~/bin_era/siesta_b4wb1preqintel18impi19fftw3 -in {self.SystemLabel}.fdf > {self.SystemLabel}.out\n")
                    # file_torque.write(f"    mpirun siesta -in {self.identifier}.fdf > {self.identifier}.out\n")
                    file_torque.write(f'    date\n')
                    file_torque.write(f'    echo "Completed fdf run"\n')
                    file_torque.write(f"    cp {self.SystemLabel}.fullBZ.WFSX {self.SystemLabel}.WFSX\n")
                    file_torque.write(f"    cp {self.SystemLabel}.selected.WFSX {self.SystemLabel}.WFSX\n")  # Since only one of them will work we dont have to worry about overwriting
                    file_torque.write(f"    ln -s ~/SIESTA_compile/siesta-master/Util/Denchar/Src/denchar .\n")
                    file_torque.write(f"    # mpirun ./denchar < {self.SystemLabel}.Denchar.fdf > {self.SystemLabel}.Denchar.out\n")
                    file_torque.write(f"    # Usage: denchar -k 1 -w 108  < *Denchar.fdf | tee Sn100_1HBDT14.denchar.out\n")
                    file_torque.write(f'    date\n')
                    # file_torque.write(f'    echo "Completed denchar run"\n')
                    # Denchar is planned to run seperatly from now on since we do not want all wavefuntions


                file_torque.write(f"END_JOB_SCRIPT\n")
            else:
                file_torque.write(f'    cd "$PWD/$dir"\n')
                file_torque.write(f"\n")
                if self.calculator == "QE":
                    file_torque.write(f"    mpirun -np {self.ntasks} {self.executable_path[self.job_handler]}pw.x -npools {self.npools} < {self.identifier}.scf.in | tee {self.identifier}.scf.out\n")
                    file_torque.write(f"    now=$(date)\n")
                    file_torque.write(f'    echo "$now : $dir : completed scf" >> ../all_jobs.log\n')
                    if "bands" in self.calculation:
                        file_torque.write(f"    mpirun -np {self.ntasks} {self.executable_path[self.job_handler]}pw.x < {self.identifier}.bands.in | tee {self.identifier}.bands.out\n")
                        file_torque.write(f"    now=$(date)\n")
                        file_torque.write(f'    echo "$now : $dir : completed bands" >> ../all_jobs.log\n')
                    if "nscf" in self.calculation:
                        file_torque.write(f"    mpirun -np {self.ntasks} {self.executable_path[self.job_handler]}pw.x < {self.identifier}.nscf.in | tee {self.identifier}.nscf.out\n")
                        file_torque.write(f"    now=$(date)\n")
                        file_torque.write(f'    echo "$now : $dir : completed nscf" >> ../all_jobs.log\n')
                    if "pdos" in self.calculation:
                        file_torque.write(f"    mpirun -np {self.ntasks} {self.executable_path[self.job_handler]}projwfc.x < {self.identifier}.pdos.in > {self.identifier}.pdos.out\n")
                        file_torque.write(f'    echo "Calculating PDOS components"\n')
                        for x in self.PDOS_required_projections:
                            if x[1] == "all":
                                file_torque.write(f"    mpirun -np {self.ntasks} {self.executable_path[self.job_handler]}sumpdos.x *\({x[0]}\)* > {self.identifier}.{x[0]}_all.PDOS\n")
                            else:
                                file_torque.write(f"    mpirun -np {self.ntasks} {self.executable_path[self.job_handler]}sumpdos.x *\({x[0]}\)*\({x[1]}\) > {self.identifier}.{x[0]}_{x[1]}.PDOS\n")
                        file_torque.write(f"    now=$(date)\n")
                        file_torque.write(f'    echo "$now : $dir : completed pdos" >> ../all_jobs.log\n')
                    if f"{self.dm}dos{self.dm}" in self.calculation:
                        file_torque.write(f'    echo "Calculationg DOS"\n')
                        file_torque.write(f"    mpirun -np {self.ntasks} {self.executable_path[self.job_handler]}dos.x < {self.identifier}.dos.in > {self.identifier}.dos.out\n")
                        file_torque.write(f"    now=$(date)\n")
                        file_torque.write(f'    echo "$now : $dir : completed dos" >> ../all_jobs.log\n')
                    if "ldos" in self.calculation:
                        file_torque.write(f'    echo "Calculationg LDOS"\n')
                        file_torque.write(f"    mpirun -np {self.ntasks} {self.executable_path[self.job_handler]}projwfc.x < {self.identifier}.ldos.in > {self.identifier}.ldos.out\n")
                        file_torque.write(f"    now=$(date)\n")
                        file_torque.write(f'    echo "$now : $dir : completed ldos" >> ../all_jobs.log\n')
                elif self.calculator == "SIESTA":
                    file_torque.write(f"    siesta -in {self.identifier}.fdf | tee {self.identifier}.out\n")
                    file_torque.write(f'    date\n')
                    file_torque.write(f'    echo "Completed fdf run"\n')
                file_torque.write("    cd .. \n")
            # file_torque.write('    cp all_jobs.log "$PWD/$dir/all_jobs.log"\n')
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
            file.write(f"#SBATCH --nodes={self.nodes}  # Number of nodes for the job\n")
            file.write(f"#SBATCH --ntasks={self.ntasks} # Number of tasks per node for the job\n")
            file.write(f"#SBATCH --mail-user=erathnayake@sivananthanlabs.us\n")
            file.write(f"#SBATCH --mail-type=ALL\n")
            if "scf" in self.calculation or "relax" in self.calculation:
                file.write(f"mpirun -np {self.ntasks} pw.x -npools {self.npools} < {self.identifier}.scf.in > {self.identifier}.scf.out\n")
                file.write(f'rm *wfc*\n')
            if "bands" in self.calculation:
                file.write(f"mpirun -np {self.ntasks} pw.x -npools {self.npools} < {self.identifier}.bands.in > {self.identifier}.bands.out\n")
            if "nscf" in self.calculation:
                file.write(f"mpirun -np {self.ntasks} pw.x -npools {self.npools} < {self.identifier}.nscf.in > {self.identifier}.nscf.out\n")
            if f"{self.dm}dos{self.dm}" in self.calculation:
                file.write(f"mpirun -np {self.ntasks} dos.x < {self.identifier}.dos.in > {self.identifier}.dos.out\n")
            if "pdos" in self.calculation:
                file.write(f"mpirun -np {self.ntasks} projwfc.x < {self.identifier}.pdos.in > {self.identifier}.pdos.out\n")
                for x in self.PDOS_required_projections:
                    if x[1] == "all":
                        file.write(f"sumpdos.x *\({x[0]}\)* > {self.identifier}.{x[0]}_all.PDOS\n")
                    else:
                        file.write(f"sumpdos.x *\({x[0]}\)*\({x[1]}\) > {self.identifier}.{x[0]}_{x[1]}.PDOS\n")
            if "ldos" in self.calculation:
                file.write(f"mpirun -np {self.ntasks} projwfc.x < {self.identifier}.ldos.in > {self.identifier}.ldos.out\n")
            # file.write(f'rm *wfc*\n')            
        os.rename(f"{self.identifier}.job", f"./{self.base_folder}/{run_name}/{self.identifier}.job")

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
                            if type(k_i) == list:
                                k_i_name = f"{k_i[0]}-{k_i[1]}-{k_i[2]}"
                            else:
                                k_i_name = f"{k_i}"
                            run_name = f"{self.identifier}{self.d}Calc{self.equals}{self.calculator}{self.d}Struct{self.equals}{self.structure_type}{self.d}Specie{self.equals}{self.specie}{self.d}KE{self.equals}{KE_cut_i}{self.d}K{self.equals}{k_i_name}{self.d}R{self.equals}{R_name}{self.d}a{self.equals}{a0_i}{self.d}type{self.equals}{self.calculation}"
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
                            if os.path.exists(f"{self.base_folder}/{run_name}"):
                                # Here we delete only folders inside the base folder. The reason is you might have multiple runs and if you delete the base folder you will delete all of them
                                shutil.rmtree(f"{self.base_folder}/{run_name}")
                                print("make_runs: Warning! Path exists!! Overwriting")
                            os.makedirs(f"{self.base_folder}/{run_name}")
                            if self.calculator == "QE":
                                self.write_QE_inputfile(run_name, KE_cut_i, R_i, a0_i, k_i)
                            elif self.calculator == "SIESTA":
                                self.write_SIESTA_inputfile(run_name)
                            self.all_runs_list.append(run_name)
                            self.move_files_to_run_folder(run_name)

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

        if sys.platform == "linux":
            extension = "sh"
            program = "#!/bin/bash\n"
        elif sys.platform == "win32":
            extension = "bat"
            program = "wsl "                

        if self.job_handler == "torque":
            bat_file = open(f"{self.base_folder}/login.{extension}", "w+")
            bat_file.write(f'{program}ssh -v -L 33301:carbon:22 rathnayake@mega.cnm.anl.gov')
            bat_file.close()
            bat_file = open(f"{self.base_folder}/session.{extension}", "w+")
            bat_file.write(f'{program}ssh -Y -p 33301 rathnayake@localhost')
            bat_file.close()
            bat_file = open(f"{self.base_folder}/rsync_out_current_folder.{extension}", "w+")
            bat_file.write(f'{program}rsync -avtuz -e "ssh -p 33301" ../{self.base_folder} rathnayake@localhost:~/Run_files/\n')
            bat_file.close()
            bat_file = open(f"{self.base_folder}/rsync_in_current_folder.{extension}", "w+")
            bat_file.write(f'{program}rsync -avtuz --exclude "*.HSX" --exclude "*.MD" --exclude "*.DM" --exclude "*.VT" --max-size=5m -e "ssh -p 33301" rathnayake@localhost:~/Run_files/{self.base_folder}/ ./\n')
            bat_file.close()
            submit_file = open(f"{self.base_folder}/run_all_jobs.sh", "w+")
            submit_file.write(f'#!/bin/bash\n')
            submit_file.write(f'dos2unix *job\n')
            submit_file.write(f'for f in *.job\n')
            submit_file.write(f'do\n')
            submit_file.write(f'  . $f\n')
            submit_file.write(f'done\n')
            submit_file.close()
            self.create_job()
        elif self.job_handler == "era_ubuntu":
            bat_file = open(f"{self.base_folder}/rsyn_out_eraubuntu.{extension}", "w+")
            bat_file.write(f'{program}rsync -avtuz -e ssh era@192.168.0.23 ./ era@192.168.0.23:~/Documents/Run_files\n')
            bat_file.close()
            bat_file = open(f"{self.base_folder}/rsyn_in_eraubuntu.{extension}", "w+")
            bat_file.write(f'{program}rsync -avtuz --max-size=5m -e ssh era@192.168.0.23:~/Documents/Run_files/ ./\n')
            bat_file.close()
            bat_file = open(f"{self.base_folder}/rsyn_in_{run_name}_eraubuntu.{extension}", "w+")
            bat_file.write(f'{program}rsync -avtuz --max-size=5m -e ssh "era@192.168.0.23:~/Documents/Run_files/{run_name}" ./\n')
            bat_file.close()
            self.create_job()
        elif self.job_handler == "slurm":
            bat_file = open(f"{self.base_folder}/rsyn_out_slurm.{extension}", "w+")
            bat_file.write(f'{program}rsync -avtuz -e ssh cluster ./ cluster:/home/erathnayake/Synced\n')
            bat_file.close()
            bat_file = open(f"{self.base_folder}/rsyn_in_slurm.{extension}", "w+")
            bat_file.write(f'{program}rsync -avtuz --max-size=5m -e ssh cluster:/home/erathnayake/Synced/ ./\n')
            bat_file.close()
            bat_file = open(f"{self.base_folder}/rsyn_in_{run_name}_slurm.{extension}", "w+")
            bat_file.write(f'{program}rsync -avtuz --max-size=5m -e ssh "cluster:/home/erathnayake/Synced/{run_name}" ./\n')
            bat_file.close()            
            # self.create_slurm_job()   already done   
        elif self.job_handler == "era_pc" or self.job_handler == "sl_laptop":
            self.create_job()  # since they are identical
        else:
            self.create_bash_file()

    def create_bash_file(self, folders = [".job.sh"], extra_pre_addlines = []):
        """
        This script creates bash files so that you can run a batch of the runs that need to be done
        """
        submit_file = open(f"{self.base_folder}/run_bash_jobs.sh", "w+")
        submit_file.write(f'#!/bin/bash\n')
        submit_file.write(f'dos2unix *job.sh\n')
        submit_file.write(f"dir_list=(")
        for x in folders:
            submit_file.write(f" {x}")
        submit_file.write(" )\n")
        # submit_file.write(f'for f in *.job.sh\n')
        submit_file.write('for f in "${dir_list[@]}"; do\n')
        submit_file.write(f'  . $f\n')
        submit_file.write(f'done\n')
        submit_file.close()
        print(f"create_bash_file: job_handler is set to: {self.job_handler}")
        bash_file = open(f"{self.base_folder}/{self.identifier}.job.sh", "w+")
        bash_file.write(f"#!/bin/bash\n\n")

        bash_file.write(f'script_start_time=$(date +%s)\n')
        bash_file.write(f'echo "Start of log" > run.log\n')
        bash_file.write(f'echo "" >> run.log\n')
        bash_file.write(f"dir_list=(")
        for x in self.all_runs_list:
            bash_file.write(f" {x}")
        bash_file.write(" )\n")
        bash_file.write(f'echo "List of folders to run" >> run.log\n')
        bash_file.write('for f in "${dir_list[@]}"; do\n')
        bash_file.write(f'  echo "$f" >> run.log\n')
        bash_file.write(f'done\n')
        bash_file.write(f'echo "" >> run.log\n\n')
        email_addresses = ",".join(self.email_recipients)
        bash_file.write(f"email_header=$'To:{email_addresses}\nFrom:statusreport_eranjan@outlook.com\nSubject:Status on: {self.identifier} Calculations\n\n'")
        bash_file.write(f"\n")
        bash_file.write(f'email_footer="\n\nOther Details\n--------------\n"\n')
        bash_file.write(f'email_footer="$email_footer Machine: $HOSTNAME\n"\n')
        bash_file.write(f'email_footer="$email_footer Solver   : {self.calculator}\n"\n')
        bash_file.write(f'email_footer="$email_footer Work Dir : $(pwd)\n\nAutomated Message\n"\n')
        bash_file.write('for dir in "${dir_list[@]}"\n')
        bash_file.write(f"do\n")
        bash_file.write(f'  cd $dir\n')
        for y in extra_pre_addlines:
            bash_file.write(f'  {y}\n')
        bash_file.write(f'  echo "Now working on: $f ... $(date)" >> ../run.log\n')
        bash_file.write(f'  dos2unix *job\n')

        if self.calculator == "SIESTA":
            bash_file.write(f'  run_start_time=$(date +%s)\n')
            bash_file.write(f'  echo "  Now working on: $f $(date)" >> ../run.log\n')
            bash_file.write(r'  mail_text="${email_header} SIESTA calculation in folder $f has started on $(date).${email_footer}"')
            bash_file.write(f"\n")
            bash_file.write(r'  echo "$mail_text" > email.txt')
            bash_file.write(f"\n")
            if self.send_mail:
                bash_file.write(f"  sendmail -t < email.txt\n")
            bash_file.write(f'  mpirun -np {self.nodes*self.procs} siesta -in {self.SystemLabel}.fdf | tee {self.SystemLabel}.out\n')
            bash_file.write(f'  run_end_time=$(date +%s)\n')
            bash_file.write(f'  elapsed_run_time=$(( run_end_time - run_start_time ))\n')
            bash_file.write(f'  mail_text="${{email_header}} Calculation in folder $f has ended on $(date). Wall_time: $elapsed_run_time s."\n')
            bash_file.write(f'  mail_text="${{email_text}} Last lines of out file:"\n')
            bash_file.write(f'  mail_text="${{email_text}} $(tail *.out)"\n')
            bash_file.write(f'  mail_text="${{email_text}}${{email_footer}}"\n')
            bash_file.write(f'  echo "$mail_text" > email_end.txt\n')
            if self.send_mail:
                bash_file.write(f'  sendmail -t < email_end.txt\n\n')
        else:
            if "scf" in self.calculation or "vc-relax" in self.calculation or "relax" in self.calculation:
                bash_file.write(f'  echo "  Now working on: scf $(date)" >> ../run.log\n')
                bash_file.write(r'  mail_text="${email_header} SCF calculation in folder $f has started on $(date).${email_footer}"')
                bash_file.write(f"\n")
                bash_file.write(r'  echo "$mail_text" > email.txt')
                bash_file.write(f"\n")
                if self.send_mail:
                    bash_file.write(f"  sendmail -t < email.txt\n")
                bash_file.write(f'  mpirun -np {self.nodes*self.ntasks} pw.x -in {self.identifier}.scf.in | tee {self.identifier}.scf.out\n')
            if "nscf" in self.calculation:
                bash_file.write(f'  echo "  Now working on: nscf $(date)" >> ../run.log\n')
                bash_file.write(r'  mail_text="${email_header} NSCF calculation in folder $f has started on $(date).${email_footer}"')
                bash_file.write(f"\n")
                bash_file.write(r'  echo "$mail_text" > email.txt')
                bash_file.write(f"\n")
                bash_file.write(f"  sendmail -t < email.txt\n")
                bash_file.write(f'  mpirun -np {self.nodes*self.ntasks} pw.x -in {self.identifier}.nscf.in | tee {self.identifier}.nscf.out\n')
            if "dos" in self.calculation:
                bash_file.write(f'  echo "  Now working on: dos $(date)" >> ../run.log\n')
                bash_file.write(r'  mail_text="${email_header} DOS calculation in folder $f has started on $(date).${email_footer}"')
                bash_file.write(f"\n")
                bash_file.write(r'  echo "$mail_text" > email.txt')
                bash_file.write(f"\n")
                bash_file.write(f"  sendmail -t < email.txt\n")
                bash_file.write(f'  mpirun -np {self.nodes*self.ntasks} pw.x -in {self.identifier}.dos.in | tee {self.identifier}.dos.out\n')
            if "pdos" in self.calculation:
                bash_file.write(f'  echo "  Now working on: pdos $(date)" >> ../run.log\n')
                bash_file.write(r'  mail_text="${email_header} PDOS calculation in folder $f has started on $(date).${email_footer}"')
                bash_file.write(f"\n")
                bash_file.write(r'  echo "$mail_text" > email.txt')
                bash_file.write(f"\n")
                bash_file.write(f"  sendmail -t < email.txt\n")
                bash_file.write(f'  mpirun -np {self.nodes*self.ntasks} pw.x -in {self.identifier}.pdos.in | tee {self.identifier}.pdos.out\n')
            if "ldos" in self.calculation:
                bash_file.write(f'  echo "  Now working on: ldos $(date)" >> ../run.log\n')
                bash_file.write(r'  mail_text="${email_header} LDOS calculation in folder $f has started on $(date).${email_footer}"')
                bash_file.write(f"\n")
                bash_file.write(r'  echo "$mail_text" > email.txt')
                bash_file.write(f"\n")
                bash_file.write(f"  sendmail -t < email.txt\n")
                bash_file.write(f'  mpirun -np {self.nodes*self.ntasks} pw.x -in {self.identifier}.ldos.in | tee {self.identifier}.ldos.out\n')
            if "bands" in self.calculation:
                bash_file.write(f'  echo "  Now working on: bands $(date)" >> ../run.log\n')
                bash_file.write(r'  mail_text="${email_header} BANDS calculation in folder $f has started on $(date).${email_footer}"')
                bash_file.write(f"\n")
                bash_file.write(r'  echo "$mail_text" > email.txt')
                bash_file.write(f"\n")
                bash_file.write(f"  sendmail -t < email.txt\n")
                bash_file.write(f'  mpirun -np {self.nodes*self.ntasks} pw.x -in {self.identifier}.bands.in | tee {self.identifier}.bands.out\n')
        bash_file.write(f'  echo "  Done" >> ../run.log\n')
        bash_file.write(r'  mail_text="${email_header} Runs in folder $f has completed on $(date).${email_footer}"')
        bash_file.write(f"\n")
        bash_file.write(r'  echo "$mail_text" > email.txt')
        bash_file.write(f"\n")
        if self.send_mail:
            bash_file.write(f"  sendmail -t < email.txt\n")
        bash_file.write(f'  cd ..\n')
        bash_file.write(f'echo "" >> run.log\n')
        bash_file.write(f'done\n')
        bash_file.write(f'\n\necho "********************************************************" >> run.log\n')
        bash_file.write(f'echo "All Runs complete. Time now is: $(date)" >> run.log\n')
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
        self.calculation      = kwargs.get("calculation", [])  # This give types of calculations in the folder eg scf+nscf+bands+dos+pdos
        self.species          = kwargs.get("species", [])
        self.high_verbosity   = kwargs.get("high_verbosity", False)
        self.MHP_base         = kwargs.get("MHP_base", "x")  # This is what you want to base you analysis on.

        # # Initializations
        self.atoms_objects = []  # Where all data read from file will be stored for a scf type file.
        self.atoms_bands_objects = []  #Where all data read from file will be stored for a bands type file.
        self.atoms_nscf_objects = []  #Where all data read from file will be stored for a nscf type file.

    def read_folder_data(self, dir):
        """
        This method reads the out files from the requried directories
        """
        # print(f"list dir method {os.listdir()}")
        # for root, dirs, files in os.walk(os.getcwd(), topdown=False):
        #     for name in dirs:
        self.directory_list = os.listdir(dir)

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

        print(f"Printing Directory list: {self.directory_list}")
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
                    if "," in x: x = x.strip("]").strip("[").split(",")  # Legacy
                    elif "-" in x: x = x.split("-")  # Current
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
            # print(self.identifier)
            # print(folder["identifier"])
            if folder["identifier"] in self.identifier or self.identifier == []:
                # print("inside Identifier")  # For debugging
                # print(f"methana")
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
                                        if folder["type"] == self.calculation or self.calculation == []:
                                            print(f"This is self.calculation:{self.calculation}")
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
                This variable sets the folder that is read to load the files.
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
        from ebk import get_machine_paths
        paths = get_machine_paths()
        if dir == "thesis":
            mydir = paths["run"]
        elif dir == "sivalabs":
            mydir = Path(runs_dir, "Run_files_SL", "Synced")
        elif dir== "here":
            mydir = Path(cur_dir)
        else:
            mydir = Path(dir)
        if self.high_verbosity:
            print(f"The Runs directory is: {mydir}")
        self.read_folder_data(mydir)
        self.make_required_folders_list()

        for x in range(0,len(self.required_folders_list)):
            path = os.path.join(mydir, self.required_folders_list[x], self.identifier[0])
            try:
                if self.folder_data[x]["Calc"].lower() == "qe":
                    # Opening Quantum Espresso Files
                    try:
                        file = ase.io.read(f"{path}.out", format = "espresso-out")  # Retains legacy code where file name might not be of the *.scf.* pattern.
                        if self.high_verbosity:
                            print(f"read_outfiles: Opening file: {path}.out")
                    except:
                        file = ase.io.read(f"{path}.scf.out", format = "espresso-out")
                        if self.high_verbosity:
                            print(f"read_outfiles: Opening file: {path}.scf.out")
                    if "bands" in self.calculation:
                        bands_file = ase.io.read(f"{path}.bands.out", format = "espresso-out")
                        if self.high_verbosity:
                            print(f"read_outfiles: Opening file: {path}.bands.out")
                    if "nscf" in self.calculation:
                        nscf_file = ase.io.read(f"{path}.nscf.out", format = "espresso-out")
                        if self.high_verbosity:
                            print(f"read_outfiles: Opening file: {path}.nscf.out")
                    if "+dos" in self.calculation:
                        self.dos_file_name = f"{path}.dos.dat"
                        if self.high_verbosity:
                            print(f"read_outfiles: Saving file name for dos file for further processing")
                elif self.folder_data[x]["Calc"].lower() == "siesta":
                    # Opening Quatnum Espresso Files
                    file = ase.io.read(f"{path}.out", format = "espresso-out")
                self.atoms_objects.append(file)
                try:
                    self.atoms_bands_objects.append(bands_file)
                except:
                    if self.high_verbosity:
                        print(f"read_outfiles: Recognized as not a bands file. bands files not appeneded to atoms_bands_objects")
                try:
                    self.atoms_nscf_objects.append(nscf_file)
                except:
                    if self.high_verbosity:
                        print(f"read_outfiles: Recognized as not a nscf file. nscf files not appeneded to atoms_nscf_objects")
            except AssertionError as er:
                print(f"read_outfiles: ** Warning Fatal Error. 'espresso.py' in ASE is giving out an assertion error as below:")
                raise
            except:
                print(f"read_outfiles: ** Warning Fatal Error. Cannot read file. File might not be present or might not have finished Recommended to set parameters to specifically exclude this file.\n{path}.out\nIt might also be the bands file of the same name.\n{path}.bands.out")
        if self.atoms_bands_objects == [] and self.atoms_nscf_objects == []:
            # No bands files have been read
            self.data = list(zip(self.required_folders_list, self.required_folder_data, self.atoms_objects))
            if self.high_verbosity:
                print(f"read_outfiles: Sucessfully read scf files. Zipping done")
        elif self.atoms_bands_objects == []:
            # Trying to zip nscf files here
            # print(self.atoms_nscf_objects)
            self.data = list(zip(self.required_folders_list, self.required_folder_data, self.atoms_objects, self.atoms_nscf_objects))
            if self.high_verbosity:
                print(f"read_outfiles: Sucessfully read nscf files and scf files. Zipping done")
        elif self.atoms_nscf_objects == []:
            # Trying to zip bands files here
            # print(self.atoms_bands_objects)
            self.data = list(zip(self.required_folders_list, self.required_folder_data, self.atoms_objects, self.atoms_bands_objects))
            if self.high_verbosity:
                print(f"read_outfiles: Sucessfully read bands files and scf files. Zipping done")
        else:
            # Trying to zip bands files here
            # print(self.atoms_bands_objects)
            self.data = list(zip(self.required_folders_list, self.required_folder_data, self.atoms_objects, self.atoms_bands_objects, self.atoms_nscf_objects))
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

    def write_structure_xyz(self, index, name):
        """
        This method writes the structure of a read file at the given 'index'.
        Index: int
            This is the index of the loaded run.
        """
        atoms = self.data[index][2]
        write(f"{name}_structure.xyz", atoms)

    def get_band_path(self):
        if self.calculation == "scf":
            print("This is a scf calculation and therefore no band path.")
        elif self.calculation == "bands":
        # try:
            return self.atoms_bands_objects[0].cell.get_band_path()
        # except:
        #     print("Error returning band path")

def make_all_job_files(base_folder = "Runs", job_list = []):
    """
    This method makes all jobs run when executing a single file. 
    Warning: Not set to properly handle .bands files since .scf has to finish in order for the .bands files to run.
             Therefore functionality is only set for .scf.job files to be listed and according to how they finish you manually run the bands.job files
    """
    print(base_folder)
    print("make_all_job_files: Printing all jobs onto a single file.")
    directory_list = os.listdir(f"{os.getcwd()}/{base_folder}")  # os.getcwd() might give different folders in different systems.
    with open(f"{base_folder}/all_jobs.job", "w+") as file:
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