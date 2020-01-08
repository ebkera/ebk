import os
import shutil
import sys
import subprocess
from ebk.kPathCreator import *
from ebk import A2Bohr

class QERunCreator:
    def __init__(self, system_name, k_set = [30,32,34], ke_set = [180, 200, 220], r_set = [300, 400]):
        """
        This function initializes the object. Default initializations are
        |self.job_manager [options: "torque", "slurm"]
        |self.nstep               : Number of molecular-dynamics or structural optimization steps performed in this run.
        |                           If set to 0, the code performs a quick "dry run", stopping just after initialization. This is useful
        |                           to check for input correctness and to have the summary printed.
        |self.tstress             : Calculate stress. It is set to .TRUE. automatically if calculation == 'vc-md' or 'vc-relax'
        |self.tprnfor             : Calculate forces. The default is false (it will be set to true automatically if calculation == 'relax','md','vc-md')
        |self.job_calculation_type: scf or bands (if bands it will do scf and bands both)
        |self.lkpoint_dir         : (Defualt: True) If .false. a subdirectory for each k_point is not opened in the "prefix".save directory; Kohn-Sham eigenvalues are
        |                           stored instead in a single file for all k-points. Currently doesn't work together with wf_collect
         
        """
        # Job initializations
        self.nodes = 1
        self.procs = 48
        self.ntasks = 48
        self.job_manager = "torque"
        self.partition = "bigmem"  # Set by default to bigmem
        self.job_calculation_type = "scf"  #options
        self.supercell_file_replace_species = []
        # QE initializations
        self.Title = "Run"
        self.system_name = system_name
        self.restart_mode = "from_scratch"
        self.k_points_number = 0
        self.k_set = k_set
        self.ke_set = ke_set
        self.r_set = r_set
        self.excorr = "pbe"     # "pbe"   = "sla+pw+pbx+pbc"    = PBE, "pz"    = "sla+pz"            = Perdew-Zunger LDA
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
        self.pseudo_dir = "./"
        self.restart_mode = "from_scratch"
        self.wf_collect = False
        self.lkpoint_dir = False
        self.disk_io = "default"
        self.verbosity  = "high"
        self.etot_conv_thr  = "1.0e-6"
        self.forc_conv_thr  = "1.0e-4"
        self.nstep  = 1
        self.tstress  = ".TRUE."
        self.tprnfor   = ".TRUE."

    def make_name(self, k, ke, r, bands, a = False):
        """
        This method creates a name for every run.
        """
        name = f"{self.system_name}_QE_K{k}_KE{ke}_R{r}"
        if a != False:
        # Here the lattice constant will only appear in the file name if you want it specifically for a lat optimization run
        # This is for another test
            name = f"{name}_a{a:.2f}"
        return name

    def jobCreator(self, k, ke, r, walltime, bands, dirname, a):
        self.walltime_days = walltime[0]
        self.walltime_hours = walltime[1]
        self.walltime_mins = walltime[2]
        self.walltime_secs = walltime[3]
        name = self.make_name(k, ke, r, bands, a)
        if bands == True:
            file_name = f"{name}.bands"
        else:
            file_name = f'{name}.scf'

        # This is if the job manager is "Torque"
        if self.job_manager == "torque":
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

        # This is if the job manager is "slurm"
        elif self.job_manager == "slurm":
            with open (f"{file_name}.job", "w") as file:
                file.write(f"#!/bin/bash\n")
                file.write(f"#SBATCH --job-name={file_name}\n")
                file.write(f"#SBATCH --partition={self.partition}\n")
                file.write(f"#SBATCH --time={self.walltime_days}-{self.walltime_hours}:{self.walltime_mins}:{self.walltime_secs}\n")
                file.write(f"#SBATCH --nodes={self.nodes}\n")
                file.write(f"#SBATCH --ntasks={self.ntasks}\n")
                file.write(f"mpirun -np {self.ntasks} pw.x -in {file_name}.in > {file_name}.out\n")
        shutil.move(f"{file_name}.job", f"./{dirname}/{file_name}.job")

    def k_file_reader(self):
        with open ("kpath.kpath",'r') as k_path_file:
            self.k_path = k_path_file.read()

    def infileCreator(self, k, ke, r, bands, dirname, a):
        # We need to have the below lines since the prefix has to not have the .scf / .bands parts
        name = self.make_name(k, ke, r, bands, a)
        # This is here to save the files names into a list for later use to call the out files
        if bands == False:
            self.set_of_files.append(name)
        if bands == True:
            file_name = f"{name}.bands"
        else:
            file_name = f'{name}.scf'
        with open (f"{file_name}.in", "w") as file:
            file.write(f"&control\n")
            file.write(f"    Title           = '{self.Title}'\n")
            file.write(f"    prefix          = '{name}'\n")
            file.write(f"    restart_mode    = '{self.restart_mode}'\n")
            file.write(f"    disk_io         = '{self.disk_io}'\n")
            if bands == False:
                file.write(f"    calculation     = 'scf'\n")
            else:
                file.write(f"    restart_mode    = '{self.restart_mode}'\n")
            file.write(f"    verbosity       = 'high'\n")
            if self.wf_collect == True:
                file.write(f"    wf_collect      = .TRUE.\n")
            if self.lkpoint_dir == False:
                file.write(f"    lkpoint_dir     = .false.\n")  # The default is true
            file.write(f"    etot_conv_thr   = {self.etot_conv_thr}\n")
            file.write(f"    forc_conv_thr   = {self.forc_conv_thr}\n")
            if self.tstress == True:
                file.write(f"    tstress         = .TRUE.\n")  # The default is false
            if self.tprnfor == True:
                file.write(f"    tprnfor         = .TRUE.\n")  # The default is false (it will be set automatically if calculation == 'relax','md','vc-md')
            if self.excorr != "auto":
                file.write(f"    input_dft       = '{self.excorr}'\n")
            if self.pseudo_dir != "system_path":
                file.write(f"    pseudo_dir      = '{self.pseudo_dir}'\n")
            file.write(f"    outdir          = './'\n")
            file.write(f"/\n")
            file.write(f"\n")
            # Start of the system block
            file.write(f"&system\n")
            file.write(f"    ibrav           = {self.ibrav}\n")
            if a != False and self.cell_parameters == False:
                # Works only if a lattice parameter has been set. and the cell_parameter_block is not set
                file.write(f"    celldm(1)       = {a}\n")
            elif a != False and self.cell_parameters == False:
                for x in range(1,6):
                    if self.celldm[x] != 0:
                        file.write(f"    celldm({x})       = {self.celldm[x]}\n")
            file.write(f"    nat             = {self.nat}\n")
            file.write(f"    ntyp            = {self.ntyp}\n")
            file.write(f"    ecutwfc         = {ke}\n")
            if bands == False:
                pass
            else:
                file.write(f"    nbnd            = 30\n")
            if self.lspinorb == True:
                file.write(f"    lspinorb        = .TRUE.\n")
                file.write(f"    noncolin        = .TRUE.\n")
            file.write(f"    ecutrho         = {r}\n")
            file.write(f"    occupations     = 'smearing'\n")
            file.write(f"    smearing        = 'gaussian'\n")
            file.write(f"    degauss         = 0.01\n")
            file.write(f"/\n")
            file.write(f"\n")
            file.write(f"&electrons\n")
            file.write(f"    mixing_beta     = 0.7\n")
            file.write(f"/\n")

            # This is for the cell_parameter_block
            if self.cell_parameters != False:
                file.write(f"\n")
                file.write(f"CELL_PARAMETERS\n")
                for vector in self.cell_parameters:
                    file.write(f"{vector[0]}   {vector[1]}   {vector[2]}\n")

            # This is for the ATOMIC_SPECIES block
            file.write(f"\n")
            file.write(f"ATOMIC_SPECIES\n")
            for specie in self.atomic_species:
                file.write(f"{specie[0]}   {specie[1]}   {specie[2]}\n")

            # This is for the ATOMIC_POSITIONS block
            file.write(f"\n")
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
                            if len(self.supercell_file_replace_species) == 0:  # No replacements set
                                if len(line.split()) == 4:
                                    line_vals = line.split()
                                    for x in range (0,len(self.supercell_file_replace_species), 2):
                                        if line_vals[0] == self.supercell_file_replace_species[x]:
                                            file.write(f"{self.supercell_file_replace_species[x+1]}    {line_vals[1]}    {line_vals[2]}    {line_vals[3]}\n")
                            else:
                                file.write(f"{line}\n")  # use this if you really dont want to do the replacement
                except:
                    print(f"QERunCreator.inFileCreator: Could not open file {self.supercell_file}. Writing default structure (diamond) with prefix Sn")
                    file.write(f"Sn        0.000000000   0.000000000   0.000000000\n")
                    file.write(f"Sn        0.250000000   0.250000000   0.250000000\n")

            if bands == False:
                file.write(f"\n")
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
            dirname = self.dirname
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
                        self.jobCreator(k, ke, r, walltime, False, dirname, b)
                        self.infileCreator(k, ke, r, False, dirname, b)
                        if self.job_calculation_type == 'bands':
                            self.infileCreator(k, ke, r, True, dirname, b)
                            self.jobCreator(k, ke, r, walltime, True, dirname, b)


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
