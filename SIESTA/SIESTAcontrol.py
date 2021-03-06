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
atomic_numbers.update({"Sn": 50})

class Generatefdf:
    def __init__(self, *args, **kwargs):
        """
        fdf_type: (string) ["bulk", "dot"]
        PAO_define: (string) ["global", "block"]-> global: Using energyshift and split norm, block: using a PAO.basis block
        PAO_define_global (bool)  # Setting this means we can set energyyshift and splitnorm. But setting this to false just means defualt value
        """
        self.SystemLabel           = kwargs.get("SystemLabel", "Sn")
        self.SystemName            = kwargs.get("SystemName", "A General Run")
        self.description           = kwargs.get("description", "A General Run")
        self.Species               = kwargs.get("Species", ["Sn", "H"])  # in order as appearing on the fdf file.
        self.coordinates_file_name = kwargs.get("coordinates_file_name", "coordinates.fdf")
        self.include_coordinate_file = kwargs.get("include_coordinate_file", False)
        self.fdf_type              = kwargs.get("fdf_type", "dot") # 
        self.PAO_define_global     = kwargs.get("PAO_define_global", False)  # Setting this so we can set energyshift and splitnorm. But setting this to false just means defualt value
        self.PAO_define            = kwargs.get("PAO_define", "global")  # SEt this to block to make the blocks work anything else and the block is ignored make sure to set the PAO_define_global to True since otherwise it will jsut be default values
        self.XC_Functional         = kwargs.get("XC_Functional", "LDA")
        self.XC_Authors            = kwargs.get("XC_Authors", "CA")
        self.LatticeVectors        = kwargs.get("LatticeVectors", [[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5,]])
        self.bands_block           = kwargs.get("bands_block", True)
        self.MPGrid                = kwargs.get("MPGrid", 10)
        self.PDOS                  = kwargs.get("PDOS", True)
        self.LDOS                  = kwargs.get("LDOS", False)
        self.PDOS_MPGrid           = kwargs.get("PDOS_MPGrid", 15)
        self.PAO_EnergyShift       = kwargs.get("PAO_EnergyShift", 0.001)
        self.PAO_SplitNorm         = kwargs.get("PAO_SplitNorm", 0.001)
        self.MD                    = kwargs.get("MD", False)
        self.Spin                  = kwargs.get("Spin", False)  # can be spin-orbit
        self.SO_strength           = kwargs.get("SO_strength", 1)
        self.include_H_in_block    = kwargs.get("include_H_in_block", False)
        self.ElectronicTemperature = kwargs.get("ElectronicTemperature", False)
        self.constrain_centre_atom = kwargs.get("constrain_centre_atom", False)  # for when doing dots we can have the central atom fixed.
        self.constrain_atom_list   = kwargs.get("constrain_atom_list", False)
        self.UseStructFile         = kwargs.get("UseStructFile", False)
        self.NetCharge             = kwargs.get("NetCharge", None)
        # Here we have all inputs for Denchar specifically
        self.Write_Denchar         = kwargs.get("Write.Denchar", False)  # This will be in the siesta fdf
        self.WriteWaveFunctions    = kwargs.get("WriteWaveFunctions", False)  # This will be in the siesta fdf
        self.Denchar_TypeOfRun     = kwargs.get("Denchar.TypeOfRun", "3D")
        self.Denchar_PlotCharge    = kwargs.get("Denchar.PlotCharge ", False)
        self.Denchar_PlotWaveFunctions = kwargs.get("Denchar.PlotWaveFunctions", True)
        self.Denchar_CoorUnits     = kwargs.get("Denchar.CoorUnits", "Ang")  # Can also be Bohr
        self.Denchar_DensityUnits  = kwargs.get("Denchar.DensityUnits ", "Ele/Ang**3")
        self.Denchar_NumberPointsX = kwargs.get("Denchar.NumberPointsX", "75")
        self.Denchar_NumberPointsY = kwargs.get("Denchar.NumberPointsX", "75")
        self.Denchar_NumberPointsZ = kwargs.get("Denchar.NumberPointsX", "75")
        self.Denchar_MinX          = kwargs.get("Denchar.MinX", "-6.5 Ang")
        self.Denchar_MaxX          = kwargs.get("Denchar.MaxX", "+6.5 Ang")
        self.Denchar_MinY          = kwargs.get("Denchar.MinY", "-6.5 Ang")
        self.Denchar_MaxY          = kwargs.get("Denchar.MaxY", "+6.5 Ang")
        self.Denchar_MinZ          = kwargs.get("Denchar.MinZ", "-6.5 Ang")
        self.Denchar_MaxZ          = kwargs.get("Denchar.MaxZ", "+6.5 Ang")  # data type is real length
        self.Denchar_PlaneGeneration = kwargs.get("Denchar.PlaneGeneration", "NormalVector")
        self.WFS_to_write_range    = kwargs.get("WFS_to_write_range", [35,40])

        if self.XC_Functional == "LDA":
            self.LatticeConstant       = kwargs.get("LatticeConstant", 6.432)
        if self.XC_Functional == "GGA":
            self.LatticeConstant       = kwargs.get("LatticeConstant", 6.672)

    def write(self, *args, **kwargs):
        with open(f"{self.SystemLabel}.fdf", "w+") as fdf_file:
            fdf_file.write(f"# -----------------------------------------------------------------------------\n")
            fdf_file.write(f"# Started on:  {datetime.now()}\n")
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
            if self.fdf_type == "bulk":
                fdf_file.write(f"LatticeConstant    {self.LatticeConstant}  Ang\t\t\t# Exchange-correlation version (PBE for GGA, PW92 or CA for LDAs)\n")

            fdf_file.write(f"\n")
            fdf_file.write(f"%block ChemicalSpeciesLabel\n")
            for i,v in enumerate(self.Species):
                fdf_file.write(f"{i+1}\t{atomic_numbers[v]}\t{v}\n")
            fdf_file.write(f"%endblock ChemicalSpeciesLabel\n")
            fdf_file.write(f"\n")
            fdf_file.write(f"MaxSCFIterations      5000\n")
            fdf_file.write(f"SCF.MustConverge      false\n")
            fdf_file.write(f"DM.MixingWeight       0.01\n\n")

            if self.fdf_type == "bulk":
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

            if self.MPGrid or self.fdf_type == "bulk":
                fdf_file.write(f"# Monkhorst-Pack Grid\n")
                fdf_file.write(f"%block kgrid.MonkhorstPack. \n")
                fdf_file.write(f"{self.MPGrid}  0  0  0.5\n")
                fdf_file.write(f"0  {self.MPGrid}  0  0.5\n")
                fdf_file.write(f"0  0  {self.MPGrid}  0.5\n")
                fdf_file.write(f"%endblock kgrid.MonkhorstPack. \n\n")

            if self.PDOS:
                if self.PDOS_MPGrid != 0:
                    fdf_file.write(f"# Monkhorst-Pack Grid\n")
                    fdf_file.write(f"%block PDOS.kgrid.MonkhorstPack. \n")
                    fdf_file.write(f"{self.PDOS_MPGrid}  0  0  0.5\n")
                    fdf_file.write(f"0  {self.PDOS_MPGrid}  0  0.5\n")
                    fdf_file.write(f"0  0  {self.PDOS_MPGrid}  0.5\n")
                    fdf_file.write(f"%endblock PDOS.kgrid.MonkhorstPack.\n\n")

                fdf_file.write(f"%block ProjectedDensityOfStates\n")
                if self.ElectronicTemperature:
                    fdf_file.write(f"-10.00 15.00 {2*self.ElectronicTemperature} 3000 eV\n")
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
                fdf_file.write(f"# These are the global values\n")
                fdf_file.write(f"PAO.BasisSize         DZP\n")
                fdf_file.write(f"PAO.EnergyShift       {self.PAO_EnergyShift} Ry   #Range of first zeta (A standard for orbital-confining cutoff radii)\n")
                fdf_file.write(f"PAO.BasisType         SPLIT      #Split Valance\n")
                fdf_file.write(f"PAO.SplitNorm         {self.PAO_SplitNorm}        #Range of second-zeta\n\n")
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

            if self.constrain_centre_atom and self.MD or self.constrain_atom_list and self.MD:
                if self.constrain_centre_atom:
                    fdf_file.write(f"%block Geometry.Constraints\n")
                    fdf_file.write(f"atom 1\n")
                    fdf_file.write(f"%endblock Geometry.Constraints\n\n")
                elif self.constrain_atom_list:
                    fdf_file.write(f"%block Geometry.Constraints\n")
                    for x in self.constrain_atom_list:
                        fdf_file.write(f"atom {x}\n")
                    fdf_file.write(f"%endblock Geometry.Constraints\n\n")

            if self.bands_block:
                fdf_file.write(f"# Band lines path\n")
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
                fdf_file.write(f"%endblock BandLines\n\n")

            if self.MD == True:
                fdf_file.write(f"MD.TypeOfRun           CG\n")
                fdf_file.write(f"MD.NumCGsteps          300\n")
                # fdf_file.write(f"MD.MaxForceTol         0.04\n")
                # fdf_file.write(f"MD.VariableCell        T\n")  # Is false by default.
                # fdf_file.write(f"MD.ConstantVolume      F\n")  # Is false by default.
                # fdf_file.write(f"MD.UseSaveXV           T\n")
                # fdf_file.write(f"MD.UseSaveCG           T\n")
                # fdf_file.write(f"MD.MaxStressTol        0.0010\n")
                # fdf_file.write(f"WriteMDHistory         T\n")
                # fdf_file.write(f"WriteMDXMol            T\n")
                # fdf_file.write(f"MD.MaxCGDispl          0.02 Bohr\n")

            fdf_file.write(f"SaveTotalPotential\t\t\t\t\ttrue\n")
            fdf_file.write(f"SaveElectrostaticPotential\ttrue\n")
            if self.UseStructFile == True:
                fdf_file.write(f"UseStructFile              true\n")
            if self.Spin: 
                fdf_file.write(f"Spin                        {self.Spin}\n")
                if self.Spin == "SpinOrbit" or self.Spin == "spin-orbit":
                    fdf_file.write(f"Spin.OrbitStrength          {self.SO_strength}\n")

            if self.NetCharge != None:
                fdf_file.write(f"NetCharge                   {self.NetCharge}\n")
            if self.ElectronicTemperature: 
                fdf_file.write(f"ElectronicTemperature       {self.ElectronicTemperature} eV\n")
            if self.Write_Denchar:
                fdf_file.write("WriteDenchar                true\n")
            if self.WriteWaveFunctions:
                fdf_file.write("COOP.Write                  true\n")
                fdf_file.write("WriteWaveFunctions          true\n")
                fdf_file.write("%block WaveFuncKPoints\n")
                fdf_file.write("0.0 0.0 0.0 from 30 to 70\n")
                fdf_file.write("%endblock WaveFuncKPoints\n")

    def write_denchar(self, *args, **kwargs):
        with open(f"{self.SystemLabel}.Denchar.fdf", "w+") as fdf_file:
            fdf_file.write(f"# -----------------------------------------------------------------------------\n")
            fdf_file.write(f"# Started on:   {datetime.now()}\n")
            fdf_file.write(f"# Description:  (Denchar) {self.description}\n")
            fdf_file.write(f"# Zincblende Sn I (alpha, grey)\n")
            fdf_file.write(f"# space group:  Fd3m\n")
            fdf_file.write(f"# lattice parameters:  a = {self.LatticeConstant} Ang ( beta tin is: a = 3.70 A, c = 3.37 A)\n")
            fdf_file.write(f"# Eranjan Kandegedara\n")
            fdf_file.write(f"# -----------------------------------------------------------------------------\n")
            fdf_file.write(f"\n")
            fdf_file.write(f"SystemName         {self.SystemName}\t\t\t\t# Descriptive name of the system\n")
            fdf_file.write(f"SystemLabel        {self.SystemLabel}\t\t\t\t# Short name for naming files\n")
            fdf_file.write(f"NumberOfSpecies    {len(self.Species)}\t\t\t\t# Number of species\n")
            fdf_file.write(f"\n")
            fdf_file.write(f"%block ChemicalSpeciesLabel\n")
            for i,v in enumerate(self.Species):
                fdf_file.write(f"{i+1}\t{atomic_numbers[v]}\t{v}\n")
            fdf_file.write(f"%endblock ChemicalSpeciesLabel\n")
            fdf_file.write(f"\n")

            fdf_file.write(f"Denchar.TypeOfRun       {self.Denchar_TypeOfRun}\n")
            if self.Denchar_PlotCharge:
                fdf_file.write(f"Denchar.PlotCharge      true\n")
            if self.Denchar_PlotWaveFunctions:
                fdf_file.write(f"Denchar.PlotWaveFunctions  true\n")

            fdf_file.write(f"Denchar.CoorUnits       {self.Denchar_CoorUnits}\n")
            fdf_file.write(f"Denchar.DensityUnits    {self.Denchar_DensityUnits}\n")
            fdf_file.write(f"\n")
            fdf_file.write(f"# Setting the mesh for Wavefunction/Charge density plot\n")
            fdf_file.write(f"Denchar.NumberPointsX   {self.Denchar_NumberPointsX}\n")
            fdf_file.write(f"Denchar.NumberPointsY   {self.Denchar_NumberPointsY}\n")
            fdf_file.write(f"Denchar.NumberPointsZ   {self.Denchar_NumberPointsZ}\n")
            fdf_file.write(f"Denchar.MinX            {self.Denchar_MinX}\n")
            fdf_file.write(f"Denchar.MaxX            {self.Denchar_MaxY}\n")
            fdf_file.write(f"Denchar.MinY            {self.Denchar_MinY}\n")
            fdf_file.write(f"Denchar.MaxY            {self.Denchar_MaxY}\n")
            fdf_file.write(f"Denchar.MinZ            {self.Denchar_MinZ}\n")
            fdf_file.write(f"Denchar.MaxZ            {self.Denchar_MaxZ}\n")


# -----------------------------------------------------------------------------
# Output options

#WriteCoorInitial     
#WriteCoorStep       
#WriteForces         
#WriteKpoints            .false.
#WriteEigenvalues        .false.
#WriteKbands             .false.
#WriteBands              .false.
#WriteMullikenPop        1            # Write Mulliken Population Analysis
#WriteCoorXmol           .false.
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