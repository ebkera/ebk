import os
import sys
import subprocess
import time
import matplotlib
import numpy as np
from datetime import datetime

atomic_numbers = {}
atomic_numbers.update({"H": 1})
atomic_numbers.update({"C": 6})
atomic_numbers.update({"N": 7})
atomic_numbers.update({"O": 8})
atomic_numbers.update({"F": 9})
atomic_numbers.update({"S": 16})
atomic_numbers.update({"Cl": 17})
atomic_numbers.update({"Ge": 32})
atomic_numbers.update({"Se": 34})
atomic_numbers.update({"Pd": 46})
atomic_numbers.update({"Sn": 50})
# Ghost atoms:
atomic_numbers.update({"H_g": -1})
atomic_numbers.update({"C_g": -6})
atomic_numbers.update({"N_g": -7})
atomic_numbers.update({"O_g": -8})
atomic_numbers.update({"F_g": -9})
atomic_numbers.update({"S_g": -16})
atomic_numbers.update({"Cl_g": -17})
atomic_numbers.update({"Ge_g": -32})
atomic_numbers.update({"Se_g": -34})
atomic_numbers.update({"Pd_g": -46})
atomic_numbers.update({"Sn_g": -50})

class Generatefdf:
    def __init__(self, *args, **kwargs):
        """
        system_type: (string) ["bulk", "dot"]  # used to be fdf_type
        PAO_define: (string) ["global", "block"]-> global: Using energyshift and split norm, block: using a PAO.basis block
        PAO_define_global (bool)  # Setting this means we can set energyyshift and splitnorm. But setting this to false just means defualt value
        """
        self.SystemLabel           = kwargs.get("SystemLabel", "Sn")
        self.SystemName            = kwargs.get("SystemName", "A General Run")
        self.description           = kwargs.get("description", "A General Run")
        self.Species               = kwargs.get("Species", ["Sn", "H"])  # in order as appearing on the fdf file.
        self.coordinates_file_name = kwargs.get("coordinates_file_name", "coordinates.fdf")
        self.include_coordinate_file = kwargs.get("include_coordinate_file", False)
        self.system_type           = kwargs.get("fdf_type", "dot") # left here for legacy code. Will be over written in teh next line if system_type is given
        self.system_type           = kwargs.get("system_type", "dot") # Can be bulk or anything else like dot or NP
        self.PAO_define_global     = kwargs.get("PAO_define_global", False)  # Setting this so we can set energyshift and splitnorm. But setting this to false just means defualt value
        self.PAO_define            = kwargs.get("PAO_define", "global")  # Set this to "block" to make the blocks work anything else and the block is ignored make sure to set the PAO_define_global to True since otherwise it will jsut be default values
        self.XC_Functional         = kwargs.get("XC_Functional", "LDA")
        self.XC_Authors            = kwargs.get("XC_Authors", "CA") # Exchange-correlation version (PBE for GGA, PW92 or CA for LDAs)
        self.LatticeVectors        = kwargs.get("LatticeVectors", [[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]])  # This only works if system_type == bulk
        self.MPGrid                = kwargs.get("MPGrid", [10,10,10])
        self.PDOS                  = kwargs.get("PDOS", True)
        self.LDOS                  = kwargs.get("LDOS", False)
        self.PDOS_MPGrid           = kwargs.get("PDOS_MPGrid", [31, 31, 31])
        self.PAO_EnergyShift       = kwargs.get("PAO_EnergyShift", None)
        self.PAO_SplitNorm         = kwargs.get("PAO_SplitNorm", None)
        self.SCF_Mix               = kwargs.get("SCF_Mix", "Hamiltonian")
        self.SCF_Mixer_Weight      = kwargs.get("SCF_Mixer_Weight", 0.01)
        self.SCF_Mixer_History     = kwargs.get("SCF_Mixer_History", 8)
        self.SCF_Mixer_Method      = kwargs.get("SCF_Mixer_Method", "Pulay")  # Options Pulay|Broyden|Linear only broyden provides steps in history
        self.SCF_Mixer_Kick        = kwargs.get("SCF_Mixer_Kick", 40)
        self.SCF_Mixer_Kick_Weight = kwargs.get("SCF_Mixer_Kick_Weight", self.SCF_Mixer_Weight)
        self.SCF_Mixer_Restart     = kwargs.get("SCF_Mixer_Restart", 0)
        self.SCF_Mixer_Restart_Save = kwargs.get("SCF_Mixer_Restart_Save", 2)
        self.SCF_Mixer_Linear_After = kwargs.get("SCF_Mixer_Linear_After", -1)  #default is -1
        self.MD                    = kwargs.get("MD", False)
        self.MD_TypeOfRun          = kwargs.get("MD_TypeOfRun", "CG")
        self.MD_Broyden_History_Steps = kwargs.get("MD_Broyden_History_Steps", 5) # Default value is 5 but set so that we remember this option
        self.MD_Steps              =kwargs.get("MD_Steps", 300)
        self.MD_MaxForceTol        =kwargs.get("MD_MaxForceTol", 0.04)  # in eV/Ang
        self.MD_VariableCell        =kwargs.get("MD_VariableCell", False)
        self.Spin                  = kwargs.get("Spin", False)  # can be spin-orbit
        self.SO_strength           = kwargs.get("SO_strength", 1)
        self.include_H_in_block    = kwargs.get("include_H_in_block", False)
        self.ElectronicTemperature = kwargs.get("ElectronicTemperature", False)
        self.constrain_centre_atom = kwargs.get("constrain_centre_atom", False)  # for when doing dots we can have the central atom fixed.
        self.constrain_atom_list   = kwargs.get("constrain_atom_list", False)
        self.constrain_cell_vectors = kwargs.get("constrain_cell_vectors", False)
        self.UseStructFile         = kwargs.get("UseStructFile", False)
        self.NetCharge             = kwargs.get("NetCharge", None)
        self.bands                 = kwargs.get("bands", True)
        self.bands_path            = kwargs.get("bands_path", None)
        self.WFS_to_write_range    = kwargs.get("WFS_to_write_range", [35,40])
        self.MaxSCFIterations = kwargs.get("MaxSCFIterations", 0)  # Here ) will set the default value (will not be written to fdf file)
        self.WriteForces           = kwargs.get("WriteForces", True)
        self.Write_Denchar         = kwargs.get("Write.Denchar", False)  # This will be in the siesta fdf
        self.WriteWaveFunctions    = kwargs.get("WriteWaveFunctions", False)  # This will be in the siesta fdf
        self.WriteMDHistory        = kwargs.get("WriteMDHistory", True)  # This will be in the siesta fdf
        self.WriteCoorStep         = kwargs.get("WriteCoorStep", False)  # This will be in the siesta fdf
        self.OpticalCalculation    = kwargs.get("OpticalCalculation", False)
        # Here we have all inputs for Denchar specifically
        self.Denchar_TypeOfRun     = kwargs.get("Denchar.TypeOfRun", "3D")
        self.Denchar_PlotCharge    = kwargs.get("Denchar.PlotCharge ", True)
        self.Denchar_PlotWaveFunctions = kwargs.get("Denchar.PlotWaveFunctions", True)
        self.Denchar_CoorUnits     = kwargs.get("Denchar.CoorUnits", "Ang")  # Can also be Bohr
        self.Denchar_DensityUnits  = kwargs.get("Denchar.DensityUnits ", "Ele/Ang**3")
        self.Denchar_NumberPointsX = kwargs.get("Denchar.NumberPointsX", "26")
        self.Denchar_NumberPointsY = kwargs.get("Denchar.NumberPointsX", "26")
        self.Denchar_NumberPointsZ = kwargs.get("Denchar.NumberPointsX", "26")
        self.Denchar_MinX          = kwargs.get("Denchar.MinX", "-6.5")
        self.Denchar_MaxX          = kwargs.get("Denchar.MaxX", "+6.5")
        self.Denchar_MinY          = kwargs.get("Denchar.MinY", "-6.5")
        self.Denchar_MaxY          = kwargs.get("Denchar.MaxY", "+6.5")
        self.Denchar_MinZ          = kwargs.get("Denchar.MinZ", "-6.5")
        self.Denchar_MaxZ          = kwargs.get("Denchar.MaxZ", "+6.5")  # data type is real length
        self.Denchar_PlaneGeneration = kwargs.get("Denchar.PlaneGeneration", "NormalVector")
        self.TS_Voltage            = kwargs.get("TS_Voltage", 0)
        self.TS_electrode_calc     = kwargs.get("TS_electrode_calculation", False)
        self.TS_ScatteringRegion_calc = kwargs.get("TS_scatteringregion_calculation", False)
        self.TS_Atoms_Buffer       = kwargs.get("TS_atoms_buffer", None)
        self.TS_ElectronicTemperature = kwargs.get("TS_ElectronicTemperature", None)
        self.TS_LelectrodePath = kwargs.get("TS_LelectrodePath", None)
        self.TS_RelectrodePath = kwargs.get("TS_RelectrodePath", self.TS_LelectrodePath)
        self.TS_LelectrodePosition = kwargs.get("TS_LelectrodePosition", 1)
        self.TS_RelectrodePosition = kwargs.get("TS_RelectrodePosition", "end -1")

        if self.XC_Functional == "LDA":
            self.LatticeConstant       = kwargs.get("LatticeConstant", 6.479)
        if self.XC_Functional == "GGA":
            self.LatticeConstant       = kwargs.get("LatticeConstant", 6.672)
        if self.XC_Functional == "VDW":
            self.LatticeConstant       = kwargs.get("LatticeConstant", 6.479)

    def write(self, *args, **kwargs):
        with open(f"{self.SystemLabel}.fdf", "w+") as fdf_file:
            fdf_file.write(f"# -----------------------------------------------------------------------------\n")
            fdf_file.write(f"# FDF file written on:  {datetime.now()}\n")
            # fdf_file.write(f"# සිංහල ජනතාව\n")
            fdf_file.write(f"# Description:  {self.description}\n")
            fdf_file.write(f"# Zincblende Sn I (alpha, grey)\n")
            fdf_file.write(f"# space group:  Fd3m\n")
            fdf_file.write(f"# lattice parameters:  a = {self.LatticeConstant} Ang ( beta tin is: a = 3.70 A, c = 3.37 A)\n")
            fdf_file.write(f"# Eranjan Kandegedara\n")
            fdf_file.write(f"# -----------------------------------------------------------------------------\n")
            fdf_file.write(f"\n")
            fdf_file.write(f"SystemName         {self.SystemName}\t\t\t\t# Descriptive name of the system\n")
            fdf_file.write(f"SystemLabel        {self.SystemLabel}\t\t\t\t# Short name for naming files\n")
            fdf_file.write(f"NumberOfSpecies    {len(self.Species)}\t\t\t\t# Number of species\n")
            fdf_file.write(f"XC.Functional      {self.XC_Functional}\t\t\t\t# Exchange-correlation functional (Defaults to LDA)\n")
            fdf_file.write(f"XC.Authors         {self.XC_Authors}\t\t\t\t# Exchange-correlation version (PBE for GGA, PW92 or CA for LDAs)\n")
            if self.system_type == "bulk":
                fdf_file.write(f"LatticeConstant    {self.LatticeConstant}  Ang\t\t\t# Exchange-correlation version (PBE for GGA, PW92 or CA for LDAs)\n")

            fdf_file.write(f"\n")
            fdf_file.write(f"%block ChemicalSpeciesLabel\n")
            for i,v in enumerate(self.Species):
                fdf_file.write(f"{i+1:<4}{atomic_numbers[v]:<4}{v:<4}\n")
            fdf_file.write(f"%endblock ChemicalSpeciesLabel\n\n")

            if self.system_type == "bulk":
                if self.LatticeVectors == "fcc":
                    fdf_file.write(f"AtomicCoordinatesFormat Fractional              # Format for coordinates\n")
                    fdf_file.write(f"%block AtomicCoordinatesAndAtomicSpecies        # Two atoms in the basis\n")
                    fdf_file.write(f".000   .000   .000   1\n")
                    fdf_file.write(f".250   .250   .250   1\n")
                    fdf_file.write(f"%endblock AtomicCoordinatesAndAtomicSpecies\n\n")
                    fdf_file.write(f"%block LatticeVectors  				#FCC lattices\n")
                    fdf_file.write(f"0.000  0.500  0.500\n")
                    fdf_file.write(f"0.500  0.000  0.500\n")
                    fdf_file.write(f"0.500  0.500  0.000\n")
                    fdf_file.write(f"%endblock LatticeVectors\n\n")
                else:
                    fdf_file.write(f"AtomicCoordinatesFormat Fractional              # Format for coordinates\n")
                    fdf_file.write(f"%block LatticeVectors  				#FCC lattices\n")
                    fdf_file.write(f"{self.LatticeVectors[0][0]}  {self.LatticeVectors[0][1]}  {self.LatticeVectors[0][2]}\n")
                    fdf_file.write(f"{self.LatticeVectors[1][0]}  {self.LatticeVectors[1][1]}  {self.LatticeVectors[1][2]}\n")
                    fdf_file.write(f"{self.LatticeVectors[2][0]}  {self.LatticeVectors[2][1]}  {self.LatticeVectors[2][2]}\n")
                    fdf_file.write(f"%endblock LatticeVectors\n\n")

            if self.include_coordinate_file:
                fdf_file.write(f"%include {self.coordinates_file_name}\n\n")
                # fdf_file.write(f"%block AtomicCoordinatesAndAtomicSpecies < {self.coordinates_file_name}\n\n")

            if self.MPGrid or self.system_type == "bulk":
                if self.MPGrid != [0,0,0]:
                    fdf_file.write(f"# Monkhorst-Pack Grid\n")
                    fdf_file.write(f"%block kgrid.MonkhorstPack. \n")
                    fdf_file.write(f"{self.MPGrid[0]}  0  0  0.5\n")
                    fdf_file.write(f"0  {self.MPGrid[1]}  0  0.5\n")
                    fdf_file.write(f"0  0  {self.MPGrid[2]}  0.5\n")
                    fdf_file.write(f"%endblock kgrid.MonkhorstPack. \n\n")

            if self.PDOS:
                if self.PDOS_MPGrid != 0:
                    fdf_file.write(f"# Monkhorst-Pack Grid - PDOS mesh\n")
                    fdf_file.write(f"%block PDOS.kgrid.MonkhorstPack. \n")
                    fdf_file.write(f"{self.PDOS_MPGrid[0]}  0  0  0.5\n")
                    fdf_file.write(f"0  {self.PDOS_MPGrid[1]}  0  0.5\n")
                    fdf_file.write(f"0  0  {self.PDOS_MPGrid[2]}  0.5\n")
                    fdf_file.write(f"%endblock PDOS.kgrid.MonkhorstPack.\n\n")

                fdf_file.write(f"%block ProjectedDensityOfStates\n")
                if self.ElectronicTemperature:
                    # fdf_file.write(f"-10.00 15.00 {2*self.ElectronicTemperature} 3000 eV\n")  # Used to be the default but now it seems that we dont get clear PDOS and more better resolution is obtained frm not taking twice the e temp
                    fdf_file.write(f"-10.00 15.00 {self.ElectronicTemperature} 3000 eV\n")
                else:
                    fdf_file.write(f"-10.00 15.00 0.050 3000 eV\n")
                fdf_file.write(f"%endblock ProjectedDensityOfStates\n\n")

            if self.LDOS:
                fdf_file.write(f"%block LocalDensityOfStates\n")
                fdf_file.write(f"-10.0  15.00 eV\n")
                fdf_file.write(f"%endblock LocalDensityOfStates\n\n")

            if self.PAO_define_global:
            # This used to be set to work only if above condition but sice we have defined the block for all we need does not hurn to have it for all other atoms
            # This way we can have H to not be forced to some value we set and have it free
                fdf_file.write(f"# Basis set optimization\n")
                fdf_file.write(f"PAO.BasisSize               DZP\n")
                if self.PAO_EnergyShift:
                    fdf_file.write(f"PAO.EnergyShift             {self.PAO_EnergyShift} Ry\t\t\t\t\t# Range of first zeta (A standard for orbital-confining cutoff radii)\n")
                fdf_file.write(f"PAO.BasisType               SPLIT    \t\t\t\t\t# Split Valance\n")
                if self.PAO_SplitNorm:
                    fdf_file.write(f"PAO.SplitNorm               {self.PAO_SplitNorm}    \t\t\t\t\t# Range of second-zeta\n\n")
            if self.PAO_define == "block":
                fdf_file.write(f"%block PAO.Basis                 # Define Basis set\n")
                # if self.XC_Functional == "LDA":
                #     fdf_file.write(f"Sn  2  # Species label, number of l-shells\n")
                #     fdf_file.write(f"  n=5  0  2  # n, l, Nzeta \n")
                #     fdf_file.write(f"  5.275  4.773\n")
                #     fdf_file.write(f"  1.000  1.000\n")
                #     fdf_file.write(f"  n=5  1  2  P  1  # n, l, Nzeta, Polarization, NzetaPol\n")
                #     fdf_file.write(f"  6.773  5.615\n")
                #     fdf_file.write(f"  1.000  1.000\n")
                # if self.XC_Functional == "GGA":
                #     fdf_file.write(f"Sn  3  # Species label, number of l-shells\n")
                #     fdf_file.write(f"  n=5  0  2  # n, l, Nzeta \n")
                #     fdf_file.write(f"  7.275,  5.143\n")
                #     fdf_file.write(f"  1.000  1.000\n")
                #     fdf_file.write(f"  n=5  1  2  # n, l, Nzeta, Polarization, NzetaPol\n")
                #     fdf_file.write(f"  7.375 4.490\n")
                #     fdf_file.write(f"  1.000  1.000\n")
                #     fdf_file.write(f"  n=5  2  1  # n, l, Nzeta, Polarization, NzetaPol\n")
                #     fdf_file.write(f"  4.170\n")
                #     fdf_file.write(f"  1.000\n")
                if "Pd" in self.Species:
                    if self.XC_Functional == "LDA":
                        fdf_file.write(f"Pd  3  # Species label, number of l-shells\n")
                        fdf_file.write(f"  n=5  0  2  # n, l, Nzeta \n")
                        fdf_file.write(f"  6.599  6.097\n")
                        fdf_file.write(f"  1.000  1.000\n")
                        fdf_file.write(f"  n=4  2  2   # n, l, Nzeta, Polarization, NzetaPol\n")
                        fdf_file.write(f"  4.735  2.517\n")
                        fdf_file.write(f"  1.000  1.000\n")
                        fdf_file.write(f"  n=5  1  1  # n, l, Nzeta, Polarization, NzetaPol\n")
                        fdf_file.write(f"  4.252\n")
                        fdf_file.write(f"  1.000\n")
                    if self.XC_Functional == "GGA":
                        fdf_file.write(f"Pd  3  # Species label, number of l-shells\n")
                        fdf_file.write(f"  n=5  0  2  # n, l, Nzeta \n")
                        fdf_file.write(f"  7.881,  3.823\n")
                        fdf_file.write(f"  1.000  1.000\n")
                        fdf_file.write(f"  n=4  2  2  # n, l, Nzeta, Polarization, NzetaPol\n")
                        fdf_file.write(f"  7.555 4.763\n")
                        fdf_file.write(f"  1.000  1.000\n")
                        fdf_file.write(f"  n=5  1  1  # n, l, Nzeta, Polarization, NzetaPol\n")
                        fdf_file.write(f"  5.474\n")
                        fdf_file.write(f"  1.000\n")
                if "Se" in self.Species:
                    if self.XC_Functional == "LDA":
                        fdf_file.write(f"Se  2  # Species label, number of l-shells\n")
                        fdf_file.write(f"  n=4  0  2  # n, l, Nzeta \n")
                        fdf_file.write(f"  4.257  3.619\n")
                        fdf_file.write(f"  1.000  1.000\n")
                        fdf_file.write(f"  n=4  1  2  P  1  # n, l, Nzeta, Polarization, NzetaPol\n")
                        fdf_file.write(f"  5.331  4.152\n")
                        fdf_file.write(f"  1.000  1.000\n")
                    if self.XC_Functional == "GGA":
                        fdf_file.write(f"Se  3  # Species label, number of l-shells\n")
                        fdf_file.write(f"  n=4  0  2  # n, l, Nzeta \n")
                        fdf_file.write(f"  7.554  4.189\n")
                        fdf_file.write(f"  1.000  1.000\n")
                        fdf_file.write(f"  n=4  1  2  # n, l, Nzeta, Polarization, NzetaPol\n")
                        fdf_file.write(f"  5.810  4.309\n")
                        fdf_file.write(f"  1.000  1.000\n")
                        fdf_file.write(f"  n=3  2  1  # n, l, Nzeta, Polarization, NzetaPol\n")
                        fdf_file.write(f"  3.708\n")
                        fdf_file.write(f"  1.000\n")
                if "C" in self.Species:
                    fdf_file.write(f"C  2  # Species label, number of l-shells\n")
                    fdf_file.write(f"  n=2  0  2  # n, l, Nzeta \n")
                    fdf_file.write(f"  6.911  3.563\n")
                    fdf_file.write(f"  1.000  1.000   \n")
                    fdf_file.write(f"  n=2  1  2  P  1  # n, l, Nzeta, Polarization, NzetaPol\n")
                    fdf_file.write(f"  9.099  3.841\n")
                    fdf_file.write(f"  1.000  1.000\n")
                if "S" in self.Species:
                    fdf_file.write(f"S  2  # Species label, number of l-shells\n")
                    fdf_file.write(f"  n=3  0  2  # n, l, Nzeta \n")
                    fdf_file.write(f"  6.702  3.587   \n")
                    fdf_file.write(f"  1.000  1.000   \n")
                    fdf_file.write(f"  n=3  1  2  P  1  # n, l, Nzeta, Polarization, NzetaPol\n")
                    fdf_file.write(f"  8.823  4.116\n")
                    fdf_file.write(f"  1.000  1.000\n")
                if "H" in self.Species:
                    fdf_file.write(f"H  1  # H from ligands with ES:0.0001\n")
                    fdf_file.write(f"  n=1  0  2  P  1  # n, l, Nzeta, Polarization, NzetaPol\n")
                    fdf_file.write(f"  8.800  4.208\n")
                    fdf_file.write(f"  1.000  1.000\n")
                fdf_file.write(f"%endblock PAO.Basis\n\n")


            fdf_file.write(f"# Electronic optimization settings\n")
            if self.Spin: 
                fdf_file.write(f"Spin                        {self.Spin}\n")
                if self.Spin == "SpinOrbit" or self.Spin == "spin-orbit":
                    fdf_file.write(f"Spin.OrbitStrength          {self.SO_strength}\n")
            if self.ElectronicTemperature: 
                fdf_file.write(f"ElectronicTemperature       {self.ElectronicTemperature} eV\n")

            fdf_file.write(f"\n# Mixing\n")
            fdf_file.write(f"SCF.Mix                     {self.SCF_Mix}\t\t\t\t\t # default: Hamiltonian others: Hamiltonian|density|charge\n")
            fdf_file.write(f"SCF.Mixer.Method            {self.SCF_Mixer_Method}\t\t\t\t\t # default: Pulay others: Pulay|Broyden|Linear\n")
            fdf_file.write(f"SCF.Mixer.Weight            {self.SCF_Mixer_Weight}\t\t\t\t\t # default: 0.25\n")
            fdf_file.write(f"SCF.Mixer.History           {self.SCF_Mixer_History}    \t\t\t\t\t # default: 2\n")
            fdf_file.write(f"SCF.Mixer.Kick              {self.SCF_Mixer_Kick}    \t\t\t\t\t # default: 0 : number of steps before linear kick to get out of local minima\n")
            fdf_file.write(f"SCF.Mixer.Kick.Weight       {self.SCF_Mixer_Kick_Weight}    \t\t\t\t\t # default: =SCF.Mixer.Weight\n")
            if self.SCF_Mixer_Restart != 0: # Removing the default case
                fdf_file.write(f"SCF.Mixer.Restart           {self.SCF_Mixer_Restart}    \t\t\t\t\t # default: =0\n")
                fdf_file.write(f"SCF.Mixer.Restart.Save      {self.SCF_Mixer_Restart_Save}    \t\t\t\t\t # default: =1 \n")
                fdf_file.write(f"SCF.Mixer.Linear.After      {self.SCF_Mixer_Linear_After}    \t\t\t\t\t # default: =-1  If≥ 0 will perform linear mixing after converged.\n")

            if self.MD == True:
                fdf_file.write(f"\n# Relaxation and Molecular Dynamics Settings\n")
                fdf_file.write(f"MD.TypeOfRun                {self.MD_TypeOfRun}\t\t\t\t\t\t # default: CG\n")
                if self.MD_TypeOfRun == "Broyden":
                    fdf_file.write(f"MD.Broyden.History.Denchar_NumberPointsZSteps    {self.MD_Broyden_History_Steps}    \t\t\t\t\t # default: 5\n")
                fdf_file.write(f"MD.Steps                    {self.MD_Steps}    \t\t\t\t # default: 0\n")
                fdf_file.write(f"MD.MaxForceTol              {self.MD_MaxForceTol} eV/Ang\t\t # default: 0.04eV/Ang\n")
                if self.MD_VariableCell:   # Requires MD.TypeOfRun to be one of the following: CG, Broyden, fire
                    fdf_file.write("MD.VariableCell             true\n")
                if self.WriteMDHistory:
                    fdf_file.write("WriteMDHistory              true\n")

                if self.constrain_centre_atom or self.constrain_atom_list or self.constrain_cell_vectors:
                    fdf_file.write(f"%block Geometry.Constraints\n")
                    if self.constrain_centre_atom:
                        fdf_file.write(f"atom 1\n")
                    elif self.constrain_atom_list:
                        for x in self.constrain_atom_list:
                            fdf_file.write(f"atom {x}\n")    # Here atom numbers start from 1
                    if self.constrain_cell_vectors:
                        for v in self.constrain_cell_vectors:
                            fdf_file.write(f"cell-vector {v}\n")    # Should be A/B/C and the vecoros have to be in the form: 100 or 010 or 001 to be considered
                    fdf_file.write(f"%endblock Geometry.Constraints\n")

            fdf_file.write(f"\n# Convergence settings\n")
            fdf_file.write(f"SCF.MustConverge            false\n")
            if self.NetCharge != None:
                fdf_file.write(f"NetCharge                   {self.NetCharge}\n")
            if self.MaxSCFIterations: 
                fdf_file.write(f"MaxSCFIterations            {self.MaxSCFIterations}\n")

            fdf_file.write(f"\n# IO settings\n")
            if self.UseStructFile == True:
                fdf_file.write(f"UseStructFile               true\n")
            fdf_file.write(f"SaveTotalPotential          true\n")
            fdf_file.write(f"SaveElectrostaticPotential  true\n")
            fdf_file.write(f"Write.DM                    true           # encouraged to have this flag true in case we have to restart or for post processing\n")
            fdf_file.write(f"MD.UseSaveXV                true\n")
            if self.WriteForces:
                fdf_file.write("WriteForces                 true\n")
            if self.Write_Denchar:
                fdf_file.write("WriteDenchar                true\n")
            if self.WriteCoorStep:
                fdf_file.write("WriteCoorStep               true\n")
            if self.WriteWaveFunctions:
                fdf_file.write("COOP.Write                  true\n")
                fdf_file.write("WriteWaveFunctions          true\n")
                fdf_file.write("%block WaveFuncKPoints\n")
                fdf_file.write("0.0 0.0 0.0 \t\t\t\t\t                     # If you want only specific EF per k point: 0.0 0.0 0.0 from 30 to 70\n")
                fdf_file.write("%endblock WaveFuncKPoints\n")

            if self.bands:
                # Use the quick bands setting to use well known band paths
                fdf_file.write(f"\n# Band structure parameters and band lines path\n")
                if self.bands_path == "2D":
                    fdf_file.write("BandLinesScale ReciprocalLatticeVectors\n")
                    fdf_file.write("%block BandLines\n")
                    fdf_file.write(" 1 0.3800000000   0.5000000000   0.0000000000     M \n")
                    fdf_file.write("40 0.0000000000   0.0000000000   0.0000000000     \GAMMA\n")
                    fdf_file.write("40 0.5000000000   0.0000000000   0.0000000000     X \n")
                    fdf_file.write("40 0.5000000000   0.5000000000   0.0000000000     S \n")
                    fdf_file.write("40 0.0000000000   0.5000000000   0.0000000000     Y \n")
                    fdf_file.write("40 0.0000000000   0.0000000000   0.0000000000     \GAMMA \n")
                    fdf_file.write("%endblock BandLines\n")

                elif self.bands_path == "FCC":
                    fdf_file.write(f"BandLinesScale pi/a\n")
                    fdf_file.write(f"%block BandLines\n")
                    fdf_file.write(f" 1  0.0000   0.0000  0.0000  \Gamma\n")
                    fdf_file.write(f"40  0.0000   2.0000  0.0000  X\n")
                    fdf_file.write(f"40  1.0000   2.0000  0.0000  W\n")
                    fdf_file.write(f"40  1.0000   1.0000  1.0000  L\n")
                    fdf_file.write(f"40  0.0000   0.0000  0.0000  \Gamma\n")
                    fdf_file.write(f"40  1.5000   1.5000  0.0000  K\n")
                    fdf_file.write(f"40  1.0000   1.0000  1.0000  L\n")
                    fdf_file.write(f"40  0.0000   0.0000  0.0000  \Gamma\n")
                    fdf_file.write(f"40  0.0000   2.0000  0.0000  X\n")
                    fdf_file.write(f"%endblock BandLines\n")

                else:
                    # Try to leave this here for default behaviour and to make a Gamma point bands file in case for reading in some values like homo and lumo edge value for a molecule calculation
                    fdf_file.write(f"BandLinesScale ReciprocalLatticeVectors\n")
                    fdf_file.write(f"%block BandLines\n")
                    fdf_file.write(f" 1  0.0000   0.0000  0.0000  \Gamma\n")
                    fdf_file.write(f"%endblock BandLines\n")

            if self.OpticalCalculation:
                # Use the quick bands setting to use well known band paths
                fdf_file.write(f"\n# Optical Property parameters\n")
                fdf_file.write(f"OpticalCalculation true\n")
                fdf_file.write(f"Optical.PolarizationType polycrystal\n")
                fdf_file.write(f"Optical.Broaden 0.2 eV\n")
                fdf_file.write("%block Optical.Mesh\n")
                fdf_file.write("  5 5 5 \n")
                fdf_file.write("%endblock Optical.Mesh\n")

            # if self.TS_electrode_calc or self.TS_ScatteringRegion_calc:
            if self.TS_electrode_calc:
                fdf_file.write(f"\n# TranSIESTA parameters\n")
                fdf_file.write(f"TS.HS.Save                 true\n")
                fdf_file.write(f"TS.DE.Save                 true\n")
            if self.TS_ScatteringRegion_calc:
                fdf_file.write(f"\n# TranSIESTA parameters\n")
                fdf_file.write(f"SolutionMethod             transiesta\n")
                # fdf_file.write(f"TS.SolutionMethod          bdt  # btd|mumps|full\n")
                fdf_file.write(f"TS.SCF.Initialize          diagon  # diagon|transiesta (avoid transiesta as still in experimental)\n")
                # fdf_file.write(f"TS.Analyze                 false\n")
                fdf_file.write("TS.Elecs.Bulk              true\n")
                fdf_file.write("TS.Elecs.DM.Update         none\n")
                fdf_file.write("TS.Elecs.GF.ReUse          false\n")
                fdf_file.write("TS.Elecs.Neglect.Principal true\n")
                # if self.TS_Voltage: fdf_file.write(f"TS.Voltage                  {self.TS_Voltage} eV\n")
                fdf_file.write(f"TS.Voltage                  {self.TS_Voltage} eV\n")
                fdf_file.write(f"TS.Hartree.Fix              +C\n")
                if self.TS_ElectronicTemperature:
                    fdf_file.write(f"TS.ElectronicTemperature   {self.TS_ElectronicTemperature}\n")
                if self.TS_Atoms_Buffer:
                    fdf_file.write(f"\n%block TS.Atoms.Buffer\n")
                    for x in self.TS_Atoms_Buffer:
                        fdf_file.write(f"atom {x}\n")    # Here atom numbers start from 1
                    fdf_file.write(f"%endblock TS.Atoms.Buffer\n")
                
                fdf_file.write(f"\n%block TS.Elecs\n")
                fdf_file.write(f"  Left\n")
                fdf_file.write(f"  Right\n")
                fdf_file.write(f"%endblock TS.Elecs\n")

                fdf_file.write(f"\n%block TS.Elec.Left\n")
                fdf_file.write(f"  HS {self.TS_LelectrodePath}\n")
                fdf_file.write(f"  chemical-potential Left\n")
                fdf_file.write(f"  semi-inf-direction -a3\n")
                fdf_file.write(f"  electrode-position {self.TS_LelectrodePosition}\n")
                fdf_file.write(f"  used-atoms 8\n")
                fdf_file.write(f"%endblock TS.Elec.Left\n")

                fdf_file.write(f"\n%block TS.Elec.Right\n")
                fdf_file.write(f"  HS {self.TS_RelectrodePath}\n")
                fdf_file.write(f"  chemical-potential Right\n")
                fdf_file.write(f"  semi-inf-direction +a3\n")
                fdf_file.write(f"  electrode-position {self.TS_RelectrodePosition}\n")
                fdf_file.write(f"  used-atoms 8\n")
                fdf_file.write(f"%endblock TS.Elec.Right\n")

                text = f"\n%block TS.ChemPots\n\
  Left\n\
  Right\n\
%endblock TS.ChemPots\n\
\n\
%block TS.ChemPot.Left\n\
  mu V/2\n\
  contour.eq\n\
    begin\n\
      c-Left\n\
      t-Left\n\
    end\n\
%endblock TS.ChemPot.Left\n\
\n\
%block TS.ChemPot.Right\n\
  mu -V/2\n\
  contour.eq\n\
    begin\n\
      c-Right\n\
      t-Right\n\
    end\n\
%endblock TS.ChemPot.Right\n\
\n\
TS.Contours.Eq.Pole 2.5 eV\n\
\n\
%block TS.Contour.c-Left\n\
  part circle\n\
   from   -40.00000 eV + V/2 to -10. kT + V/2\n\
    points 30\n\
     method g-legendre\n\
%endblock TS.Contour.c-Left\n\
\n\
%block TS.Contour.t-Left\n\
  part tail\n\
   from prev to inf\n\
    points 10\n\
     method g-fermi\n\
%endblock TS.Contour.t-Left\n\
\n\
%block TS.Contour.c-Right\n\
  part circle\n\
   from   -40.00000 eV - V/2 to -10. kT - V/2\n\
    points 30\n\
     method g-legendre\n\
%endblock TS.Contour.c-Right\n\
\n\
%block TS.Contour.t-Right\n\
  part tail\n\
   from prev to inf\n\
    points 10\n\
     method g-fermi\n\
%endblock TS.Contour.t-Right\n\
\n\
TS.Contours.nEq.Eta 0.0001 eV\n\
\n\
%block TS.Contours.nEq\n\
  neq\n\
%endblock TS.Contours.nEq\n\
\n\
%block TS.Contour.nEq.neq\n\
  part line\n\
   from -|V|/2 - 5 kT to |V|/2 + 5 kT\n\
    delta 0.01 eV\n\
     method mid-rule\n\
%endblock TS.Contour.nEq.neq\n\
"
                fdf_file.write(text)
                # fdf_file.write(f"MD.MaxForceTol         0.04\n")
                # fdf_file.write(f"MD.VariableCell        T\n")  # Is false by default.
                # fdf_file.write(f"MD.ConstantVolume      F\n")  # Is false by default.
                # fdf_file.write(f"MD.UseSaveXV           T\n")
                # fdf_file.write(f"MD.UseSaveCG           T\n")
                # fdf_file.write(f"MD.MaxStressTol        0.0010\n")
                # fdf_file.write(f"WriteMDHistory         T\n")
                # fdf_file.write(f"WriteMDXMol            T\n")
                # fdf_file.write(f"MD.MaxCGDispl          0.02 Bohr\n")            

    # fdf_file.write(f"\n# Denchar settings\n")
    def write_denchar(self, *args, **kwargs):
        with open(f"{self.SystemLabel}.Denchar.fdf", "w+") as fdf_file:
            fdf_file.write(f"# -----------------------------------------------------------------------------\n")
            fdf_file.write(f"# Started on:   {datetime.now()}\n")
            fdf_file.write(f"# Description:  (Denchar) {self.description}\n")
            fdf_file.write(f"# Zincblende Sn I (alpha, grey)\n")
            fdf_file.write(f"# space group:  Fd3m\n")
            fdf_file.write(f"# lattice parameters:  a = {self.LatticeConstant} Ang ( beta tin is: a = 3.70 A, c = 3.37 A)\n")
            fdf_file.write(f"# Eranjan Kandegedara\n\n")
            fdf_file.write(f"# Usage: denchar -k 1 -w 108  < *Denchar.fdf | tee Sn100_1HBDT14.denchar.out\n")
            fdf_file.write(f"#        Where -k : kpoint of interest, -w : wavefunction of interest\n")
            fdf_file.write(f"# -----------------------------------------------------------------------------\n")
            fdf_file.write(f"\n")
            fdf_file.write(f"SystemName         {self.SystemName:<15}\t\t\t\t# Descriptive name of the system\n")
            fdf_file.write(f"SystemLabel        {self.SystemLabel:<15}\t\t\t\t# Short name for naming files\n")
            fdf_file.write(f"NumberOfSpecies    {len(self.Species):<15}\t\t\t\t# Number of species\n")
            fdf_file.write(f"\n")
            fdf_file.write(f"%block ChemicalSpeciesLabel\n")
            for i,v in enumerate(self.Species):
                fdf_file.write(f"{i+1:<4}\t{atomic_numbers[v]:<4}\t{v:<4}\n")
            fdf_file.write(f"%endblock ChemicalSpeciesLabel\n")
            fdf_file.write(f"\n")

            fdf_file.write(f"Denchar.TypeOfRun           {self.Denchar_TypeOfRun:<15}\n")
            if self.Denchar_PlotCharge:
                fdf_file.write(f"Denchar.PlotCharge          true\n")
            if self.Denchar_PlotWaveFunctions:
                fdf_file.write(f"Denchar.PlotWaveFunctions   false\n")

            fdf_file.write(f"Denchar.CoorUnits           {self.Denchar_CoorUnits:<15}\n")
            fdf_file.write(f"Denchar.DensityUnits        {self.Denchar_DensityUnits:<15}\n")
            fdf_file.write(f"\n")
            fdf_file.write(f"# Setting the mesh for Wavefunction/Charge density plot\n")
            fdf_file.write(f"Denchar.NumberPointsX       {self.Denchar_NumberPointsX:<15}\n")
            fdf_file.write(f"Denchar.NumberPointsY       {self.Denchar_NumberPointsY:<15}\n")
            fdf_file.write(f"Denchar.NumberPointsZ       {self.Denchar_NumberPointsZ:<15}\n")
            fdf_file.write(f"Denchar.MinX                {self.Denchar_MinX:<15} Ang\n")
            fdf_file.write(f"Denchar.MaxX                {self.Denchar_MaxY:<15} Ang\n")
            fdf_file.write(f"Denchar.MinY                {self.Denchar_MinY:<15} Ang\n")
            fdf_file.write(f"Denchar.MaxY                {self.Denchar_MaxY:<15} Ang\n")
            fdf_file.write(f"Denchar.MinZ                {self.Denchar_MinZ:<15} Ang\n")
            fdf_file.write(f"Denchar.MaxZ                {self.Denchar_MaxZ:<15} Ang\n")


# -----------------------------------------------------------------------------
# Output options

#WriteCoorInitial     
#WriteForces         
#WriteKpoints            .false.
#WriteEigenvalues        .false.
#WriteKbands             .false.
#WriteBands              .false.
#WriteMullikenPop        1            # Write Mulliken Population Analysis
#WriteMDCoorXmol         .false.
#WriteMDhistory          .false.
#WriteCoorXmol           .false.

# -----------------------------------------------------------------------------
# Options for saving/reading information

#DM.UseSaveDM                         # Use DM Continuation files
#MD.UseSaveXV               .false.   # Use stored positions and velocities
#MD.UseSaveCG               .false.   # Use stored positions and velocities
#SaveRho                              # Write valence pseudocharge at the mesh
#SaveDeltaRho                         # Write RHOscf-RHOatm at the mesh
#SaveElectrostaticPotential .false.   # Write the total elect. pot. at the mesh
#SaveTotalPotential         .false.   # Write the total pot. at the mesh
#WriteSiestaDim             .false.   # Write minimum dim to siesta.h and stop
#WriteDenchar                         # Write information for DENCHAR


# -----------------------------------------------------------------------------
# Notes 
#
# Things Chi-Gong wanted me to check for.
# The valance orbitals - mine are 4f10 5s2 5p2 5d0 as per the paper we are following
# The cutoff radii for the particular orbitals and we are using values from the paper
# The The number of Zetas.. We are using double zeta and also split valance.
# Karanna thiyenne, 
# Get the bond lengths vs Energy optimization
# Get the cohesive energies.
#     - Cohesive energy gives us the energy per bond basically
# The the double derivative of the min E vs a0.
#     - This gives us the bulk modulus

if __name__ == "__main__":
    pass

    # This is an example and is commented out so that nothing actually happens if this file is accidentally run and "pass" is set above this to ensure that.
    # # making the fdf_file
    # kwargs = {"SystemLabel":label, "fdf_type" : "NP_ligand", "PAO_define" : "block", "ElectronicTemperature" : 0.0244, "MD": True, "PDOS": True, "PDOS_MPGrid": 0, "LDOS": False, "include_coordinate_file": True}
    # if XC == "LDA": kwargs.update({"XC_Functional":"LDA", "XC_Authors":"CA"})
    # kwargs.update({"NumberOfSpecies": 4})
    # kwargs.update({"coordinates_file_name": f"{label}_coordinates.fdf"})
    # kwargs.update({"constrain_centre_atom" : False})
    # kwargs.update({"constrain_atom_list": atoms_to_constrain})
    # kwargs.update({"Spin": "spin-orbit"})
    # kwargs.update({"MPGrid" : 0})
    # kwargs.update({"SO_strength": SO})
    # kwargs.update({"bands_block": False})
    # kwargs.update({"PAO_define_global": True})
    # kwargs.update({"PAO_SplitNorm": 0.30})
    # kwargs.update({"include_H_in_block": True})