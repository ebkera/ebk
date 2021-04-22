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
    shutil.copy(f"{RELAX_DIR}/CONTCAR", f"{folder_name}/POSCAR")

    
