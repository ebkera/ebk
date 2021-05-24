import os
import shutil


def make_folder_structure():
    try:
        os.mkdir("1_RELAX")
        os.mkdir("2_SCF")
        os.mkdir("3_DOS")
        os.mkdir("4_BANDS")
        os.mkdir("DEFAULTS")
    except OSError as error:
        print(error)

def populate_RELAX(folder_name = False, RELAX_DIR="1_RELAX"):
    if folder_name == False:
        folder_name = RELAX_DIR
    else:
        folder_name = f"{RELAX_DIR}_{folder_name}"
    
    if not os.path.isdir(folder_name):
        os.mkdir(folder_name)
    shutil.copy(f"DEFAULTS/POTCAR", f"{folder_name}/POTCAR")
    shutil.copy(f"DEFAULTS/POSCAR_unrelaxed", f"{folder_name}/POSCAR")
    # shutil.copy(f"DEFAULTS/INCAR__TEMPLATE", f"{folder_name}/INCAR")
    # Lets make the INCAR file
    file = open(f"{folder_name}/INCAR", "w+")
    contents = get_relaxation_INCAR()
    file.write(contents)
    # Lets make the KPOINTS file
    file = open(f"{folder_name}/KPOINTS", "w+")
    contents = get_relaxation_KPOINTS()
    file.write(contents)
    file.close()

    # shutil.copy(f"DEFAULTS/KPOINTS_SCF_TEMPLATE", f"{folder_name}/KPOINTS")
    # shutil.copy(f"{RELAX_DIR}/CONTCAR", f"{folder_name}/POSCAR")

def populate_SCF(folder_name = False, RELAX_DIR="1_RELAX"):
    if folder_name == False:
        folder_name = "2_SCF"
    else:
        folder_name = f"2_SCF_{folder_name}"
    
    if not os.path.isdir(folder_name):
        os.mkdir(folder_name)
    shutil.copy(f"DEFAULTS/POTCAR", f"{folder_name}/POTCAR")
    try:
        shutil.copy(f"DEFAULTS/INCAR_SCF_TEMPLATE", f"{folder_name}/INCAR")
        print(f"INCAR for SCF run loaded from template file in DEFAULTS")
    except:
        print(f"No template for INCAR for SCF using original template from code")
        file = open(f"{folder_name}/INCAR", "w+")
        contents = get_scf_INCAR()
        file.write(contents)
        file.close()
    try:
        shutil.copy(f"DEFAULTS/KPOINTS_SCF_TEMPLATE", f"{folder_name}/KPOINTS")
        print(f"KPOINTS for SCF run loaded from template file in DEFAULTS")
    except:
        print(f"No template for KPOINTS for SCF using original template from code")
        file = open(f"{folder_name}/KPOINTS", "w+")
        contents = get_scf_KPOINTS()
        file.write(contents)
        file.close()
    try:
        shutil.copy(f"{RELAX_DIR}/CONTCAR", f"{folder_name}/POSCAR")
        print(f"CONTCAR file found in relaxation run.. Using it for POSCAR in SCF")
    except:
        print(f"No CONTCAR file found in relaxation run.. Looking for DEFAULTS/POSCAR_relaxed.. if you already have relaxed structure rename as such..")
        shutil.copy(f"DEFAULTS/POSCAR_relaxed", f"{folder_name}/POSCAR")

def populate_BANDS(folder_name = False, RELAX_DIR="1_RELAX", SCF_DIR="2_SCF"):
    if folder_name == False:
        folder_name = "4_BANDS"
    else:
        folder_name = f"4_BANDS_{folder_name}"
    
    os.mkdir(folder_name)
    shutil.copy(f"DEFAULTS/POTCAR", f"{folder_name}/POTCAR")
    shutil.copy(f"DEFAULTS/INCAR_BANDS_TEMPLATE", f"{folder_name}/INCAR")
    shutil.copy(f"DEFAULTS/KPOINTS_BANDS_TEMPLATE", f"{folder_name}/KPOINTS")
    # shutil.copy(f"{SCF_DIR}/CHGCAR", f"{folder_name}/CHGCAR")
    shutil.copy(f"{RELAX_DIR}/CONTCAR", f"{folder_name}/POSCAR")

def populate_epsilon(folder_name = False, RELAX_DIR="1_RELAX", SCF_DIR="2_SCF"):
    folder_base = "4_eps2"
    if folder_name == False:
        folder_name = folder_base
    else:
        folder_name = f"{folder_base}_{folder_name}"
    
    os.mkdir(folder_name)
    shutil.copy(f"DEFAULTS/POTCAR", f"{folder_name}/POTCAR")
    shutil.copy(f"DEFAULTS/INCAR_eps_TEMPLATE", f"{folder_name}/INCAR")
    shutil.copy(f"DEFAULTS/KPOINTS_DOS_TEMPLATE", f"{folder_name}/KPOINTS")
    # shutil.copy(f"{SCF_DIR}/CHGCAR", f"{folder_name}/CHGCAR")
    shutil.copy(f"{RELAX_DIR}/CONTCAR", f"{folder_name}/POSCAR")


def make_NSCF_calculation(folder_list, run_name="run", SCF_DIR="2_SCF", out_file =f"run.log"):
    folder_list_text = ''
    for x in folder_list:
        folder_list_text+=(f' "{x}"')

    string_to_write = f'#!/bin/bash\n\
\n\
folder_list=({folder_list_text} )\n\
\n\
echo "Start of log" > {out_file}\n\
echo "" >> {out_file}\n\
echo "List of folders to run" >> {out_file}\n\
\n\
for f in "${{folder_list[@]}}"; do\n\
    echo "$f" >> {out_file}\n\
done\n\
\n\
echo "" >> {out_file}\n\
\n\
for f in "${{folder_list[@]}}"; do\n\
    cd $f\n\
    echo "Now working on $f ... $(date)" >> ../{out_file}\n\
    # cp ../2_SCF/KPOINTS\n\
    cp ../{SCF_DIR}/CHGCAR CHGCAR\n\
    cp ../{SCF_DIR}/POSCAR POSCAR\n\
    cp ../{SCF_DIR}/POTCAR POTCAR\n\
    mpirun -np 4 vasp_std_NON_SO | tee era.out\n\
    cd ..\n\
done\n\
echo "done"\n\
echo "done" >> {out_file}'

    with open(f"{run_name}.sh", "w+") as file:
        file.write(string_to_write)


def get_relaxation_INCAR():
    content = f"SYSTEM = RELAXATION_for_\n\
  \n\
# start parameters for this Run (automatic defaults are finem, hence not often required)\n\
  ISTART = 1         # job   : 0-new  1-orbitals from WAVECAR (continuation job:restart with constant energy cut-off) 3-orbitals from WAVECAR (continuation job: restart with constant basis set)\n\
  ICHARG = 2         # charge: 1: from CHGCAR file | 2-atom (for SCF) | 10+: NSCF calculations\n\
  PREC   = Accurate  # standard precision (OtherOptions: Accurate)\n\
  \n\
# electronic optimization\n\
  ENCUT  = 420.00 eV # defaults from POTCAR, but wise to include (Recommened 130% of POTCAT val)\n\
  ALGO   = Normal    # alorithm for electron optimization, can be also FAST or ALL\n\
  NELM   = 200       # of ELM steps, sometimes default is too small \n\
  EDIFF  = 1E-07     # stopping-criterion for ELM\n\
  ISMEAR = 0         # -5:Tetrahedral Method for smearing \n\
  SIGMA  = 0.05\n\
  # ENMAX  = 400       # cutoff should be set manually  (This seems to be an obsolete flag...)\n\
  # AMIN   = 0.01      # Default: 0.10 specifies the minimal mixing parameter in Kerker's[1] initial approximation to the charge dielectric function used in the Broyden[2][3]/Pulay[4] mixing scheme (IMIX=4, INIMIX=1)\n\
  # LSORBIT = .TRUE.   # Spin Orbit Coupling is set to true.\\\n\
  \n\
  # van der Waals\n\
  # IVDW    = 1         # IVDW=1|10 DFT-D2 method of Grimme (available as of VASP.5.2.11)\n\
  # VDW_RADIUS=50.0     # cutoff radius (in Å {{\displaystyle \AA }} \AA ) for pair interactions\n\
  # VDW_S6  =0.75 	    # global scaling factor s 6 {{\displaystyle s_{{6}}}} s_{{6}} (available in VASP.5.3.4 and later)\n\
  # VDW_SR  =		    # 1.00 scaling factor s R {{\displaystyle s_{{R}}}} s_{{R}} (available in VASP.5.3.4 and later)\n\
  # VDW_D   =20.0 	    # damping parameter d {{\displaystyle d}} d\n\
  # VDW_C6=[real array] C 6 \n\
  # VDW_R0=[real array] R 0 \n\
  # LVDW_EWALD=.FALSE.  # decides whether lattice summation in E d i s p\n\
  \n\
# DOS calculations\n\
  LORBIT = 11        # 11 for both total and projected\n\
  NEDOS  = 1000      # numbr of points for DOS\n\
  #  EMIN   = -5\n\
  #  EMAX   = 5\n\
  #  NBANDS = *\n\
  \n\
# ionic relaxation\n\
  ISIF   = 0         # 0: only atoms nostress | 1: Relaxing atoms stress trace only | 2: Relaxing atoms stress trace full | 4:cell shape, and cell volume\n\
  EDIFFG = -1E-02    # stopping-criterion for IOM (If negative: all forces smaller 1E-2)\n\
  NSW    = 200       # number of steps for IOM in other words 20 ionic steps\n\
  POTIM  = .5        # step for ionic-motion (for MD in fs)\n\
  NFREE  = 2         # depending on IBRION, NFREE specifies the number of remembered steps in the history of ionic convergence runs, or the number of ionic displacements in frozen phonon calculations. however systems of low dimensionality require a careful setting of NFREE (or preferably an exact counting of the number of degrees of freedom)\n\
  IBRION = 1         # -1:no update. | 0:molecular dynamics.| 1:ionic relaxation (RMM-DIIS) (usually faster) | 2:ionic relaxation (conjugate gradient algorithm) | 3 | 5 | 6 | 7 | 8 | 44  \n\
  \n\
# performance optimization\n\
  NCORE   = 2         # one orbital handled by 4 cores recommened: 4-(number of cores)\n\
  #  LREAL  = A        # real space projection; slightly less accurate but faster \n\
  #  KPAR   = 2        # make 4 groups, each group working on one set of k-points "
	
    # print(content)
    return content

def get_relaxation_KPOINTS():
    content = f"Automatic Mesh\n\
 0\n\
Monkhorst Pack\n\
 5 5 1\n\
 0. 0. 0.\n\
 "
    return content

def get_scf_KPOINTS():
    content = f"Automatic Mesh\n\
 0\n\
Monkhorst Pack\n\
 15 15 1\n\
 0. 0. 0.\n\
 "
    return content

def get_scf_INCAR():
    content = f"SYSTEM = SCF_for_\n\n\
# start parameters for this Run (automatic defaults are finem, hence not often required)\n\
  ISTART = 0         # job   : 0-new  1- orbitals from WAVECAR\n\
  ICHARG = 2         # charge: 1: from CHGCAR file | 2-atom (for SCF) | 10+: NSCF calculations\n\
  PREC   = Accurate  # standard precision (OtherOptions: Accurate)\n\
\n\
# electronic optimization\n\
  ENCUT  = 420.00 eV # defaults from POTCAR, but wise to include (Recommened 130% of POTCAT val)\n\
  ALGO   = Normal    # alorithm for electron optimization, can be also FAST or ALL\n\
  NELM   = 200       # of ELM steps, sometimes default is too small \n\
  EDIFF  = 1E-07     # stopping-criterion for ELM\n\
  ISMEAR = 0         # -5:Tetrahedral Method for smearing \n\
  SIGMA  = 0.05\n\
  # ENMAX  = 400       # cutoff should be set manually  (This seems to be an obsolete flag...)\n\
  # AMIN   = 0.01      # Default: 0.10 specifies the minimal mixing parameter in Kerker's[1] initial approximation to the charge dielectric function used in the Broyden[2][3]/Pulay[4] mixing scheme (IMIX=4, INIMIX=1)\n\
  # LSORBIT = .TRUE.   # Spin Orbit Coupling is set to true.\\\n\
  \n\
  # van der Waals\n\
  # IVDW    = 1         # IVDW=1|10 DFT-D2 method of Grimme (available as of VASP.5.2.11)\n\
  # VDW_RADIUS=50.0     # cutoff radius (in Å {{\displaystyle \AA }} \AA ) for pair interactions\n\
  # VDW_S6  =0.75 	    # global scaling factor s 6 {{\displaystyle s_{{6}}}} s_{{6}} (available in VASP.5.3.4 and later)\n\
  # VDW_SR  =		    # 1.00 scaling factor s R {{\displaystyle s_{{R}}}} s_{{R}} (available in VASP.5.3.4 and later)\n\
  # VDW_D   =20.0 	    # damping parameter d {{\displaystyle d}} d\n\
  # VDW_C6=[real array] C 6 \n\
  # VDW_R0=[real array] R 0 \n\
  # LVDW_EWALD=.FALSE.  # decides whether lattice summation in E d i s p\n\
\n\
# DOS calculations\n\
  LORBIT = 11        # 11 for bith total and projected\n\
  NEDOS  = 1000      # numbr of points for DOS\n\
#  EMIN   = -5\n\
#  EMAX   = 5\n\
#  NBANDS = *\n\
\n\
# ionic relaxation\n\
  ISIF   = 0         # 0: only atoms nostress | 1: Relaxing atoms stress trace only | 2: Relaxing atoms stress trace full | 4:cell shape, and cell volume\n\
  EDIFFG = -1E-02    # stopping-criterion for IOM (If negative: all forces smaller 1E-2)\n\
  NSW    = 0         # number of steps for IOM in other words 20 ionic steps\n\
  POTIM  = .5        # step for ionic-motion (for MD in fs)\n\
  NFREE  = 2         # 2 independent degrees of freedom\n\
  IBRION = -1         # IBRION = -1:no update. | 0:molecular dynamics.| 1:ionic relaxation (RMM-DIIS) (usually faster) | 2:ionic relaxation (conjugate gradient algorithm) | 3 | 5 | 6 | 7 | 8 | 44  \n\
\n\
# performance optimization\n\
  NCORE   = 2         # one orbital handled by 4 cores recommened: 4-(number of cores)\n\
#  LREAL  = A        # real space projection; slightly less accurate but faster \n\
#  KPAR   = 2        # make 4 groups, each group working on one set of k-points "
	
    # print(content)
    return content