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

def populate_SCF(folder_name = False, RELAX_DIR="1_RELAX"):
    if folder_name == False:
        folder_name = "2_SCF"
    else:
        folder_name = f"2_SCF_{folder_name}"
    
    os.mkdir(folder_name)
    shutil.copy(f"DEFAULTS/POTCAR", f"{folder_name}/POTCAR")
    shutil.copy(f"DEFAULTS/INCAR_SCF_TEMPLATE", f"{folder_name}/INCAR")
    shutil.copy(f"DEFAULTS/KPOINTS_SCF_TEMPLATE", f"{folder_name}/KPOINTS")
    shutil.copy(f"{RELAX_DIR}/CONTCAR", f"{folder_name}/POSCAR")

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

