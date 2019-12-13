from ebk.kPathCreator import *
from ebk.QEBandsReader import QEBandsReader
from ebk.BandPlotter import BandPlotter
import os
# import matplotlib.pyplot as plt
# matplotlib.use('Agg')  # no UI backend required if working in the wsl without a UI (Un comment this if on WSL)
import sys
import subprocess

def run(*args, **kwargs):
    #Defaults are set here
    name = 'run'
    tee = True
    bands = False
    np = 1
    for key, value in kwargs.items():
        if key == 'name': name = value
        if key == 'tee': tee = value
        if key == 'bands': bands = value
        if key == 'np': np = value

    if tee == True:
        out_type = "| tee"
    else:
        out_type = ">"

    if bands == True:
        bands = ".bands"
    else:
        bands = ".scf"

    print(f"EBK: {bands} run with {np} processors")
    os.system(f"mpirun  -np {np} '/mnt/c/Users/Eranjan/Desktop/Quantum_Expresso/qe-6.4.1/bin/pw.x' -in {name}{bands}.in {out_type} {name}{bands}.out")

if __name__ == "__main__":
    run (name = "Sn", np = 3, bands = False, tee = True)
    # run (name = "Sn", np = 3, bands = True, tee = True)

    path = kPathCreator()
    path.add_startval(G)
    path.add_kpath(X, 20)
    path.add_kpath(W, 20)
    path.add_kpath(L, 20)
    path.add_kpath(G, 20)
    path.add_kpath(K, 20)
    path.add_kpath(L, 20)

    # path.add_startval(L)
    # path.add_kpath(G, 30)
    # path.add_kpath(X, 30)
    # path.add_kpath(U, 30)
    # path.add_kpath(G, 30)
    # path.out_kpath_QE()

    # myplot = QEBandsReader("Sn.bands.out")
    # plotter= BandPlotter(path.k_distance, myplot.data, path.highSymPoints_position, path.highSymPoints_symbols)
    # plotter.title = "Sn QE Band Diagram"
    # plotter.file_name = "Sn_K20_KE63_R300"
    # plotter.vlines = True
    # plotter.set_y_range = True
    # plotter.ylim_high = 5
    # plotter.ylim_low = -10
    # plotter.same_band_colour = False
    # plotter.Ef_shift =  8.9383
    # plotter.plot()

    # plotter2= BandPlotter(path.k_distance, myplot.data, path.highSymPoints_position, path.highSymPoints_symbols)
    # plotter2.title = "Sn QE Band Diagram"
    # # plotter2.file_name = "Ef9shift"
    # plotter2.vlines = True
    # plotter2.set_y_range = False
    # plotter2.same_band_colour = True
    # plotter2.Ef_shift = 0
    # plotter2.band_colour = "r"
    # # plotter2.new_fig = True
    # plotter2.plot()