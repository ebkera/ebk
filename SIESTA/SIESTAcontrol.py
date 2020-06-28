import os
import sys
import subprocess
import time
import matplotlib
import numpy as np
from datetime import datetime

class Generatefdf:
    def __init__(self, *args, **kwargs):
        """
        fdf_type: (string) ["bulk", "dot"]
        PAO_define: (string) ["global", "block"]-> global: Using energyshift and split norm, block: using a PAO.basis block
        """
        self.SystemLabel           = kwargs.get("SystemLabel", "Sn")
        self.lattice_contant       = kwargs.get("lattice_contant", 6.432)
        self.description           = kwargs.get("description", "A General Run")
        self.NumberOfSpecies       = kwargs.get("NumberOfSpecies", 2)
        self.Species               = kwargs.get("Species", ["Sn", "H"])
        self.coordinates_file_name = kwargs.get("coordinates_file_name", "coordinates.fdf")
        self.fdf_type              = kwargs.get("fdf_type", "dot") # 
        self.PAO_define            = kwargs.get("PAO_define", "global")
        self.XC_Functional         = kwargs.get("XC.Functional", "LDA")
        self.XC_Authors            = kwargs.get("XC", "CA")
        self.lattice_vectors       = kwargs.get("lattice_vectors", "fcc")

    def write(self, *args, **kwargs):
        with open(f"{self.SystemLabel}.fdf", "w+") as fdf_file:
            fdf_file.write(f"# -----------------------------------------------------------------------------\n")
            fdf_file.write(f"# Started on:  {datetime.now()}\n")
            fdf_file.write(f"# Description:  {self.description}\n")
            fdf_file.write(f"# Zincblende Sn I (alpha, grey)\n")
            fdf_file.write(f"# space group:  Fd3m\n")
            fdf_file.write(f"# lattice parameters:  a = {self.lattice_contant} A ( beta tin is: a = 3.70 A, c = 3.37 A)\n")
            fdf_file.write(f"# Eranjan Kandegedara\n")
            fdf_file.write(f"# -----------------------------------------------------------------------------\n")
            fdf_file.write(f"\n")
            fdf_file.write(f"\n")
            fdf_file.write(f"SystemName         alpha-Sn\t\t\t\t# Descriptive name of the system\n")
            fdf_file.write(f"SystemLabel        {self.SystemLabel}\t\t\t\t# Short name for naming files\n")
            if self.fdf_type == "bulk":
                fdf_file.write(f"NumberOfSpecies    1\t\t\t\t# Number of species\n")
            elif self.fdf_type == "dot":
                fdf_file.write(f"NumberOfSpecies    2\t\t\t\t# Number of species\n")
            else self.fdf_type == "other":
                fdf_file.write(f"NumberOfSpecies    {self.NumberOfSpecies}\t\t\t\t# Number of species\n")
            fdf_file.write(f"XC.Functional      {self.XC_Functional}\t\t\t\t# Exchange-correlation functional (Defaults to LDA)\n")
            fdf_file.write(f"XC.Authors         {self.XC_Authors}\t\t\t\t# Exchange-correlation version (PBE for GGA, PW92 or CA for LDAs)\n")
            fdf_file.write(f"\n")
            fdf_file.write(f"%block Chemical_Species_Label\n")
            fdf_file.write(f"1   50    Sn\n")
            if self.fdf_type == "dot": fdf_file.write(f"2    1    H\n")
            # fdf_file.write(f"3    6    C\n")
            # fdf_file.write(f"4   16    S\n")
            fdf_file.write(f"%endblock Chemical_Species_Label\n")
            fdf_file.write(f"\n")
            fdf_file.write(f"MaxSCFIterations      5000\n")
            fdf_file.write(f"SCF.MustConverge      false\n")
            fdf_file.write(f"DM.MixingWeight       0.01\n")
            fdf_file.write(f"\n")
            if self.fdf_type == "bulk":
                if self.lattice_vectors == "fcc":
                    fdf_file.write(f"%block LatticeVectors  				#FCC lattices\n")
                    fdf_file.write(f"0.000  0.500  0.500\n")
                    fdf_file.write(f"0.500  0.000  0.500\n")
                    fdf_file.write(f"0.500  0.500  0.000\n")
                    fdf_file.write(f"%endblock LatticeVectors\n\n")
            if self.fdf_type != "bulk":
                fdf_file.write(f"%include {self.coordinates_file_name}\n\n")
            fdf_file.write(f"%block ProjectedDensityOfStates\n")
            fdf_file.write(f"-10.00 15.00 00.05 3000 eV\n")
            fdf_file.write(f"%endblock ProjectedDensityOfStates\n\n")
            fdf_file.write(f"# These values are from the paper\n")
            fdf_file.write(f"PAO.BasisSize         DZP\n")
            fdf_file.write(f"PAO.EnergyShift       0.001 Ry    #Range of first zeta (A standard for orbital-confining cutoff radii)\n")
            fdf_file.write(f"PAO.BasisType         SPLIT       #Split Valance\n")
            fdf_file.write(f"PAO.SplitNorm         0.30        #Range of second-zeta\n")
            fdf_file.write(f"\n")
            if self.PAO_define == "block":
                fdf_file.write(f"%block PAO.Basis                 # Define Basis set\n")
                fdf_file.write(f"Sn  2  # Species label, number of l-shells\n")
                fdf_file.write(f"  n=5  0  2  # n, l, Nzeta \n")
                fdf_file.write(f"  5.275  4.773\n")
                fdf_file.write(f"  1.000  1.000\n")
                fdf_file.write(f"  n=5  1  2  P  1  # n, l, Nzeta, Polarization, NzetaPol\n")
                fdf_file.write(f"  6.773  5.615\n")
                fdf_file.write(f"  1.000  1.000\n")
                # fdf_file.write(f"C  2  # Species label, number of l-shells\n")
                # fdf_file.write(f"  n=2  0  2  # n, l, Nzeta \n")
                # fdf_file.write(f"  6.911  3.563\n")
                # fdf_file.write(f"  1.000  1.000   \n")
                # fdf_file.write(f"  n=2  1  2  P  1  # n, l, Nzeta, Polarization, NzetaPol\n")
                # fdf_file.write(f"  9.099  3.841\n")
                # fdf_file.write(f"  1.000  1.000\n")
                # fdf_file.write(f"S  2  # Species label, number of l-shells\n")
                # fdf_file.write(f"  n=3  0  2  # n, l, Nzeta \n")
                # fdf_file.write(f"  6.702  3.587   \n")
                # fdf_file.write(f"  1.000  1.000   \n")
                # fdf_file.write(f"  n=3  1  2  P  1  # n, l, Nzeta, Polarization, NzetaPol\n")
                # fdf_file.write(f"  8.823  4.116\n")
                # fdf_file.write(f"  1.000  1.000\n")
                if self.fdf_type == "dot":
                    fdf_file.write(f"H  1  # H from Sn dots.\n")
                    fdf_file.write(f"  n=1  0  2  P  1  # n, l, Nzeta, Polarization, NzetaPol\n")
                    fdf_file.write(f"  7.026  3.359\n")
                    fdf_file.write(f"  1.000  1.000\n")
                fdf_file.write(f"%endblock PAO.Basis\n")
            fdf_file.write(f"\n")
            fdf_file.write(f"SaveTotalPotential true\n")
            fdf_file.write(f"SaveElectrostaticPotential true\n")


if __name__ == "__main__":
    run = optimize_parameter()
