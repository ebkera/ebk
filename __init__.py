"""
This program has all the constants and device parameters that we will need and also some 
conversion functions that is very commonly used in code in this package.
It will also include functions that we would have to use all the time but does not belong 
anywhere else.
More diverse funtionality has been added since initial planned scope was set. 
"""

#Here are some global innitializations
# import os
import sys
from numpy import ALLOW_THREADS, positive
import scipy.constants
import math

# from matplotlib import pyplot as plt
import matplotlib
font = {'family' : "Times New Roman",
        # 'weight' : 'bold',
        # 'size'   : 22
        }
matplotlib.rc('font', **font)
# plt.rcParams["font.family"] = "Times New Roman"
# Setting matplotlib to load non gui back end if not on a windows platform THis shoudl only be enabled for wsl systems
# print(sys.platform)
if sys.platform != "win32" and sys.platform != "linux":
    print("ebk: System might not have a GUI backend. Suspecting system is WSL on Windows")
    matplotlib.use('Agg')  # no UI backend required if working in the wsl without a UI

# Universal constants
Egap = 0.0
c = 2.99792458E8
e = scipy.constants.e # 1.60218E-19
m0 = scipy.constants.electron_mass # 9.109E-31  # electron mass in Kg
h = 6.62607015E-34  # (Units Js) Starting from 2019 May 20 This value has been fixed
heV = 4.135667662E-15  # (Units eVs) This is the value is eVs
epsilon0 = scipy.constants.epsilon_0 #8.85E-12 # (Units F/m) 

# Bulk masses of material here
me_CdSe = 1.18E-31
mh_CdSe = 4.09E-31
me_HgTe = 4E-4*m0
mh_HgTe = 0.5*m0
me_Sn = 0.0236*m0  # alpha-Sn
mh_Sn = 0.21*m0  # alpha-Sn

# Lattice Constants here in Angstroms
a_InAs = 6.0583

def get_machine_paths():
    import os
    """This function will get the machine paths to the directories that contain teh following depending on the computer in use.
    1) Pseudopotentials
    2) structure files (eg. xyz)
    3) Path to the run files (run)
    This function returns a dict with the relevant paths
    """
    # print(os.uname()[1])  # for testing
    try:
        if os.uname()[1] == "ERA-PC":
            print("get_machine_paths: loading path for ERA-PC")
            pps = f"C:/Users/Eranjan/OneDrive - University of Illinois at Chicago/CQD Research/Run_files/PseudopotentialDatabase"
            xyz = f"C:/Users/Eranjan/OneDrive - University of Illinois at Chicago/CQD Research/Run_files/XYZdatabase"
            run = f"C:/Users/Eranjan/OneDrive - University of Illinois at Chicago/CQD Research/Run_files_git/Run_files"
        if os.uname()[1] == "era":
            print("get_machine_paths: loading path for era-ubuntu")
            pps = f"/home/era/OneDrive_UIC/CQD Research/Run_files/PseudopotentialDatabase"
            xyz = f"/home/era/OneDrive_UIC/CQD Research/Run_files/XYZdatabase"
            run = f"/home/era/OneDrive_UIC/CQD Research/Run_files_git/Run_files"
        if os.uname()[1] == "ffl4":
            print("get_machine_paths: loading path for ffl4: probably Pop! OS")
            pps = f"/home/ebk/OneDrive_UIC/CQD Research/Run_files/PseudopotentialDatabase"
            xyz = f"/home/ebk/OneDrive_UIC/CQD Research/Run_files/XYZdatabase"
            run = f"/home/ebk/OneDrive_UIC/CQD Research/Run_files_git/Run_files"
    except:
            print(f"Defaulting to SivaLabs variables")
            pps = f"/mnt/c/Users/erathnayake/OneDrive - University of Illinois at Chicago/CQD Research/Run_files/PseudopotentialDatabase"
            xyz = f"/mnt/c/Users/erathnayake/OneDrive - University of Illinois at Chicago/CQD Research/Run_files/XYZdatabase"
            run = f"/mnt/c/Users/erathnayake/OneDrive - University of Illinois at Chicago/CQD Research/Run_files_git/Run_files"
    paths = {"pps": pps, "xyz": xyz, "run": run}
    return paths

def reduced_mass(me2, mh2):
    '''This function calculates reduced mass'''
    return me2*mh2/(me2+mh2)

def J2mum(E_J):
    return 1e6*h*c/E_J

def J2nm(E_J):
    return 1e9*h*c/E_J

def eV2mum(E_eV):
    return 1e6*h*c/eV2J(E_eV)

def eV2nm(E_eV):
    return 1e9*h*c/eV2J(E_eV)

def mum2eV(mum):
    return J2eV(h*c/(mum*1E-6))

def nm2eV(mum):
    return J2eV(h*c/(mum*1E-9))

def J2eV(J):
    return J/e

def eV2J(eV):
    return eV*e

def A2Bohr(A):
    """Takes distances in Angstroms and returns distances in Bohr"""
    return A*1.88973

def Bohr2A(B):
    """Takes distances in Bohr and returns distances in Angstroms"""
    return B*0.529177

def Rydberg2eV(R):
    """Takes energies in Rydbergs and returns energies in eV"""
    return R*13.6056980659

def eV2Rydberg(e):
    """Takes energies in eV and returns energies in Rydbergs"""
    return e/13.6056980659

def Rydberg2J(R):
    """Takes energies in Rydbergs and returns energies in J"""
    return R*2.179872E-18

def J2Rydberg(J):
    """Takes energies in J and returns energies in Rydbergs"""
    return J*4.5874253167158E+17 

def eVA32GPa(e):
    """Takes eV/Angstroms^3 and converts it into GPa. This is useful when calculating the bulk modulus especially."""
    return e*160.21766208

def gap_calculator_1D(m_e, m_h, L):
    return J2eV((scipy.constants.pi**2*scipy.constants.hbar**2)/(2*m_e*(L*1E-9)**2) + (scipy.constants.pi**2*scipy.constants.hbar**2)/(2*m_h*(L*1E-9)**2))

def eV2K(eV):
    """Takes in electronic temperature in eV and outputs in K"""
    return eV*1.160E4

def K2eV(K):
    """Takes in electronic temperature in K and outputs in eV"""
    return K/1.160E4

def get_pseudopotential(id):
    """
    This method returns the pseudopotential that corresponds to the identifier that I have set so that we can quickly identify a psedo. This is just for internal tracking
    Obviously these identifiers will be only required for a pseudopotential optimization run.
    """
    identifier = id
    print(identifier)
    if identifier == "ONCV_PBE_FR1.1":
        pseudopotential = {'Sn': f"Sn_ONCV_PBE_FR-1.1.upf"}
    elif identifier == "ONCV_PBE_SR1.0":
        pseudopotential = {'Sn': f"Sn_ONCV_PBE-1.0.upf"}
    elif identifier == "ONCV_PBE_SR1.1":
        pseudopotential = {'Sn': f"Sn_ONCV_PBE-1.1.upf"}
    elif identifier == "ONCV_PBE_SR1.2":
        pseudopotential = {'Sn': f"Sn_ONCV_PBE-1.2.upf"}
    elif identifier == "Sn_slab_FR":
        pseudopotential = {'Sn': f"Sn_ONCV_PBE_FR-1.1.upf", "H":"H_ONCV_PBE_FR-1.0.upf"}
    elif identifier == "Sn_slab_SR1.1":
        pseudopotential = {'Sn': f"Sn_ONCV_PBE-1.1.upf", "H":"H_ONCV_PBE_FR-1.0.upf"}
    elif identifier == "Sn_slab_SR1.2":
        pseudopotential = {'Sn': f"Sn_ONCV_PBE-1.2.upf", "H":"H_ONCV_PBE_FR-1.0.upf"}
    elif identifier == "Sn_slabs":
        pseudopotential = {'Sn': f"Sn_ONCV_PBE_FR-1.1.upf", "H":"H_ONCV_PBE_FR-1.0.upf"}
    # PBE pps
    # FR
    elif identifier == "QE_PBE_FR_NLCC_1":
        pseudopotential = {'Sn': f"Sn.rel-pbe-dn-kjpaw_psl.1.0.0.UPF"}
    elif identifier == "QE_PBE_FR_NLCC_2":
        pseudopotential = {'Sn': f"Sn.rel-pbe-dn-rrkjus_psl.1.0.0.UPF"}
    elif identifier == "QE_PBE_FR_NLCC_3":
        pseudopotential = {'Sn': f"Sn.rel-pbe-dn-kjpaw_psl.0.2.UPF"}
    elif identifier == "QE_PBESOL_FR_NLCC_1":
        pseudopotential = {'Sn': f"Sn.rel-pbesol-dn-kjpaw_psl.1.0.0.UPF"}
    # PBE pps
    # SR
    elif identifier == "QE_PBE_SR_NLCC_1":
        pseudopotential = {'Sn': f"Sn.pbe-dn-kjpaw_psl.1.0.0.UPF"}
    elif identifier == "QE_PBE_SR_NLCC_2":
        pseudopotential = {'Sn': f"Sn.pbe-dn-rrkjus_psl.1.0.0.UPF"}
    elif identifier == "QE_PBE_SR_NLCC_3":
        pseudopotential = {'Sn': f"Sn.pbe-dn-kjpaw_psl.0.2.UPF"}
    elif identifier == "QE_PBESOL_SR_NLCC_1":
        pseudopotential = {'Sn': f"Sn.pbesol-dn-kjpaw_psl.1.0.0.UPF"}
    #LDA pps
    elif identifier == "QE_PZ_FR_NLCC_1":
        pseudopotential = {'Sn': f"Sn.rel-pz-dn-kjpaw_psl.0.2.UPF"}
    elif identifier == "MCT_Slab_EDT" or identifier == "MCT_Slab_BDT12" or identifier == "MCT_Slab_BDT14":
        pseudopotential = {'C': f"c_pbe_v1.2.uspp.F.UPF", "H":"h_pbe_v1.4.uspp.F.UPF", "Hg":"hg_pbe_v1.uspp.F.UPF", "S":"s_pbe_v1.4.uspp.F.UPF", "Te":"te_pbe_v1.uspp.F.UPF"}
    elif "Slabs_passivation" in identifier:
        # This is for the relaxation of slabs
        pseudopotential = {'C': f"c_pbe_v1.2.uspp.F.UPF", "H":"h_pbe_v1.4.uspp.F.UPF", "Hg":"hg_pbe_v1.uspp.F.UPF", "Te":"te_pbe_v1.uspp.F.UPF"}
    elif identifier == "Slabs_passivation_TeTop" or identifier == "Slabs_passivation_HgTop" or identifier == "MCT_Slab":
        # This is for the relaxation of slabs
        pseudopotential = {"H":"h_pbe_v1.4.uspp.F.UPF", "Hg":"hg_pbe_v1.uspp.F.UPF", "Te":"te_pbe_v1.uspp.F.UPF"}
        #For the ligand: pps are ultrasoft
    elif identifier == "EDT" or identifier == "BDT12" or identifier == "BDT14":
        pseudopotential = {'C': 'c_pbe_v1.2.uspp.F.UPF', 'H': 'h_pbe_v1.4.uspp.F.UPF', 'S': 's_pbe_v1.4.uspp.F.UPF'}
    elif "MCT_Slab_110_bilayers_" in identifier:
        # This is for the relaxation of slabs
        pseudopotential = {"H":"h_pbe_v1.4.uspp.F.UPF", "Hg":"hg_pbe_v1.uspp.F.UPF", "Te":"te_pbe_v1.uspp.F.UPF"}         
    elif "MCT_Slab_with_ligands" in identifier:
        # This is for the relaxation of slabs
        pseudopotential = {"H":"h_pbe_v1.4.uspp.F.UPF", "Hg":"hg_pbe_v1.uspp.F.UPF", "Te":"te_pbe_v1.uspp.F.UPF", 'C': 'c_pbe_v1.2.uspp.F.UPF', 'H': 'h_pbe_v1.4.uspp.F.UPF', 'S': 's_pbe_v1.4.uspp.F.UPF'}                 
    return pseudopotential

def make_mini_proj():
    """
    Useful when wanting to create a quick mini project
    """
    import os
    os.mkdir("tex_media")
    tex_file_text = r"""%This is a mini project
\documentclass[a4paper,12pt]{article}
\usepackage[hidelinks]{hyperref}
\usepackage{float}
\usepackage{graphicx}
\graphicspath{{tex_media/}} %Setting the graphicspath

\title{Sn Bulk and NP band/energy structure optimization - with tuned DZP basis set}
\author{Eranjan Kandegedara}
\begin{document}
\maketitle
\tableofcontents

\section{Background}
This is dummy text.
\begin{figure}[H]
    \centering
    \includegraphics[scale = 0.45]{Rivero_band_gaps.JPG}
    \caption{Bands from the rivero paper. We are looking at the dashed green line.}
    \label{fig:rivero_bands}
\end{figure}

\bibliographystyle{unsrt}
\bibliography{../../referencesCQD}
\end{document}"""

    with open("main.tex", "w+") as tex_file:
        tex_file.write(tex_file_text)
        
def add_xyz_files(file1, file2, file_name_to_write=f"combined.xyz"):
    """
    Adds two .xyz files.
    Input: two .xyz files without the extension
    Output: Single .xyz file with all the atoms in the same file.
            See code for naming scheme for similar atoms.
    """

    # Saving input files to variables
    file = open(f"{file1}.xyz", 'r')
    data_file1 = [line for line in file]
    file.close()
    file = open(f"{file2}.xyz", 'r')
    data_file2 = [line for line in file]
    file.close()

    # The number of atoms here.
    N_file1 = int(data_file1[0])
    N_file2 = int(data_file2[0])
    N_total = N_file1 + N_file2

    file = open(f"{file_name_to_write}", 'w+')
    file.write(f"{N_total}\n")
    file.write("This is a combined file for comparison of two structures - Eranjan\n")

    for x in range(2,len(data_file1)):
        line = data_file1[x].split()
        file.write(f"A{line[0]}\t{float(line[1]):.7f}\t{float(line[2]):.7f}\t{float(line[3]):.7f}\n")

    for x in range(2,len(data_file2)):
        line = data_file2[x].split()
        file.write(f"B{line[0]}\t{float(line[1]):.7f}\t{float(line[2]):.7f}\t{float(line[3]):.7f}\n")
    file.close()    

def subtract_xyz_files(file1, file2, file_name_to_write=f"combined.xyz"):
    """
    Subtracts two .xyz files.
    Input: two .xyz files without the extension
    Output: Single .xyz file with all the atoms in the same file.
            See code for naming scheme for similar atoms.
    """

    # Saving file contents to variables
    file = open(f"{file1}.xyz", 'r')
    data_file1 = [line for line in file]
    file.close()
    file = open(f"{file2}.xyz", 'r')
    data_file2 = [line for line in file]
    file.close()

    # The number of atoms here.
    N_file1 = int(data_file1[0])
    N_file2 = int(data_file2[0])

    f1_lines = []
    f2_lines = []
    
    # Here we convert the lines into a list and also makes the numerical values into floats so we can compare 
    for x in range(2,len(data_file1)):
        line = data_file1[x].split()
        f1_lines.append([f"A{line[0]}", float(line[1]), float(line[2]), float(line[3])])

    for x in range(2,len(data_file2)):
        line = data_file2[x].split()
        f2_lines.append([f"B{line[0]}", float(line[1]), float(line[2]), float(line[3])])

    final_f1_lines = []
    final_f2_lines = []

    for line in f1_lines:
        save = True
        for l in f2_lines:
            if abs(line[1] - l[1]) <= 0.001 and abs(line[2] - l[2]) <= 0.001 and abs(line[3] - l[3]) <= 0.001:
                save = False
                # print(line, l)
                break
        if save:
            final_f1_lines.append(line)
        
    for line in f2_lines:
        save = True
        for l in f1_lines:
            if abs(line[1] - l[1]) <= 0.001 and abs(line[2] - l[2]) <= 0.001 and abs(line[3] - l[3]) <= 0.001:
                save = False
                # print(line, l)
                break
        if save:
            final_f2_lines.append(line)

    N_total = len(final_f1_lines) + len(final_f2_lines)

    file = open(f"{file_name_to_write}", 'w+')
    file.write(f"{N_total}\n")
    file.write("This is a combined (subtracted) file for comparison of two structures - Eranjan\n")

    for x in range(0,len(final_f1_lines)):
        line = final_f1_lines[x]
        # print(line)
        file.write(f"{line[0]}\t{float(line[1]):.7f}\t{float(line[2]):.7f}\t{float(line[3]):.7f}\n")

    for x in range(0,len(final_f2_lines)):
        line = final_f2_lines[x]
        file.write(f"{line[0]}\t{float(line[1]):.7f}\t{float(line[2]):.7f}\t{float(line[3]):.7f}\n")
    file.close() 

def make_ghost(file_list, atoms_list, output_file_name="out.xyz"):
    """
    Inputs:
        file_list: (list) a list of files were atoms in atoms_list will be made into ghost atoms by renaming the species. Without extension.
        atoms_list: (list of lists) : a list with length equal to file_list, where each element is a list of atoms to ghost. Can also be keywords like "all"
    This will change chemical symbols all atoms in the atoms_list.
    """
    file_contents_list = []
    for file in file_list:
        file_ = open(f"{file}.xyz", 'r')
        file_contents_list.append([line for line in file_])
        file_.close()

    # print(file_contents_list)
    total_number_of_atoms = 0
    for file_contents in file_contents_list:
        total_number_of_atoms += int(file_contents[0])

    file = open(f"{output_file_name}", 'w+')
    file.write(f"{total_number_of_atoms}\n")
    file.write("This file has been modified to contain ghost atoms labeled ChemSym_g- Eranjan\n")
    for i,file_contents in enumerate(file_contents_list):

        # for x in range(0,len(file_contents_list)):
        for j,line in enumerate(file_contents):
            # print(line)
            if j <= 1: 
                print("continuing")
                continue
            line = line.split()
            # print(atoms_list[i][j])
            if atoms_list[i] == "all" or "all" in atoms_list[i] or (j-2) in atoms_list[i]:
                file.write(f"{line[0]}_g\t{float(line[1]):.7f}\t{float(line[2]):.7f}\t{float(line[3]):.7f}\n")
            else:
                file.write(f"{line[0]}\t{float(line[1]):.7f}\t{float(line[2]):.7f}\t{float(line[3]):.7f}\n")


def get_species(atoms):
    """
    returns the species as a list for use in other scripts
    input should be an atoms type object 
    maybe we need to implement for an xyz type file as well
    """

    from ase.atoms import Atoms as atoms_type
    species=[]

    if type(atoms) is atoms_type:
        all_symbols = atoms.get_chemical_symbols()
        for x in all_symbols:
            if x not in species:
                species.append(x)

    return species


class progress_bar():
    """To be used as a progress bar"""    
    def __init__(self, tasks, descriptor = ""):
        self.resolution = 100
        self.tasks = tasks
        self.d = tasks/self.resolution
        self.descriptor = descriptor
        # print(self.d)

    def get_progress(self, completed):
        current_bucket = 0
        # We usually expect inputs to be indeces list or otherwise but the number of tasks to be the length of the list so we increment by 1
        completed+=1
        # if (completed % self.d) == 0:

        prog = math.ceil(completed/self.d)
        if prog == current_bucket: return
        else:current_bucket = prog
        # print(prog)
        to_prog = self.resolution-prog

        p=f" {self.descriptor}|"
        # p=f" |"
        for x in range(0,prog-1):
            p+="="
        p+=">"
        # p+=f"{prog}%"
        for x in range(to_prog):
            p+=" "
        p+="|"
        
        p = f"{p}{prog}% ({completed}/{self.tasks})"

        if (completed/self.d) == 100 or completed >= self.tasks:
            print(f"{p}")
        else:
            print(f"{p} ", end="\r")

