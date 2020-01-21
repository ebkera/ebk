import os
import sys
import subprocess
import time
import matplotlib
matplotlib.use('Agg')  # no UI backend required if working in the wsl without a UI
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np

class Optimization:
    def __init__(self):
        print(__name__)

    class OptimizationLatticeConstant(Optimization):
        def __init__(self):
            




# startTime = time.time()	
# # Set parameters here
# a_0 = [6.20,6.22,6.24,6.26,6.28,6.30,6.32,6.34,6.36,6.38,6.40,6.432,6.457,6.475,6.50,6.55,6.60,6.672]
# #a_0 = [6.20,6.22]
# output = "a" # a or v
# dirname = 'a0OptimizationRun_MP10(.5)_2'

# # These are the energy values
# E_vals = []

# # Making at directory to put all the diagrams and out files in.
# os.mkdir(dirname)

# for x in a_0:
#     # Get all the lines from the current fdf file.
#     file = open("Sn.fdf", "r")
#     rows = []
#     for line in file:
#         rows.append(line)
#     file.close
    
#     # Save new fdf file
#     for row in rows:
#         if "LatticeConstant" in row:
#             index = rows.index(row)
#             rows[index] = "LatticeConstant        "+str(x)+"  Ang     # 6.457 in paper for no rel corrs and 6.432 for final values and 6.672 for GGA final\n"
    
#     # Overwriting the fdf file
#     file_toWrite = open("Sn.fdf", "w+")
#     for row in rows:
#         file_toWrite.write(row)
#     file_toWrite.close()
    
#     # Run plot file which will run siesta
#     os.system("python3 plot.py r g h 0 1 2 3 4 5 6 7 8 9 10")

#     # Read the Sn.out file
#     file = open("Sn.out", "r")
#     rows = []
#     for line in file:
#         rows.append(line)
#     file.close()
#     # saving the Energy values in a list
#     for row in rows:
#         if "siesta:         Total =" in row:
#             E_vals.append(row.strip("siesta:         Total ="))
    
#     #Waiting a bit for the plots to be saved into a pdf just to be be safe
#     print("Waiting a bit for the plots to be saved into a pdf just to be be safe")
#     time.sleep(5)
#     #Renaming the band diagrams
#     os.rename("Sn_BandDiagramPlots.pdf", "./"+dirname+"/Sn_BandDiagramPlots_" + str(x) +".pdf")	
#     #Renaming the outfiles
#     os.rename("Sn.out", "./"+dirname+"/Sn" + str(x) +".out")
#     #Renaming the bands files
#     os.rename("Sn.bands", "./"+dirname+"/Sn" + str(x) +".bands")
    

# ## Calculations and plotting the lattice constants
# to_plot = []
# for i in E_vals:
#     to_plot.append(float(i))

# # Since we will need volume/ per atom
# v_0 = [n**3/4.0 for n in a_0]

# # calculate polynomial
# z = np.polyfit(a_0, to_plot, 3) #abnits
# v = np.polyfit(v_0, to_plot, 3) #abnits

# #z = np.polyfit(a_01, E_LDA_RELwkgrid_Monkhorst_Pack6, 3) #relativistic
# f = np.poly1d(z)
# f_v = np.poly1d(v)

# # calculate new x's and y's
# x_new = np.linspace(a_0[0], a_0[-1], 100)
# Vx_new = np.linspace(v_0[0], v_0[-1], 100)
# y_new = f(x_new)
# Vy_new = f_v(Vx_new)

# # Getting the minimum energy
# E_optimized = min(y_new)
# E_optimized_printable = E_optimized.astype(np.float)

# # Getting the optimized lattice constant
# min_index = np.argmin(y_new)
# min_index_v = np.argmin(Vy_new)
# a0_optimized = x_new.flat[min_index]
# v0_optimized = Vx_new.flat[min_index_v]

# # Write Energies to an energies.my file
# file_toWrite_output = open("./"+dirname+"/energies.my", "w+")
# for x in E_vals:
#     file_toWrite_output.write(x)		
# file_toWrite_output.write("The optimized lattice constant is: " + str(a0_optimized) + "\n")
# file_toWrite_output.write("The optimized volume is: " + str(v0_optimized) + "\n")
# file_toWrite_output.write("The optimized Total Energy is: " + str(E_optimized_printable) + "\n")
# file_toWrite_output.write("The fit parameters: \n" )
# #z_printable = z.astype(np.float)
# #file_toWrite_output.write(z_printable)
# file_toWrite_output.close()

# # Calculations
# # Getting the double derivative using a 3rd degree polynomial
# dda0 = 6*z.flat[0]*a0_optimized + 2*z.flat[1]
# ddv0 = 6*v.flat[0]*v0_optimized + 2*v.flat[1]
# B = v0_optimized*ddv0*160.21766208 # 1 eV/Angstrom3 = 160.21766208 GPa

# # The Individual plots.
# if output == "a":
#     toplotx = a_0
#     fitx = x_new
#     fity = y_new
#     plt.xlabel('Lattice Constant a$_0$ in $\\AA$')
# elif output == "v":
#     toplotx = v_0
#     fitx = Vx_new
#     fity = Vy_new
#     plt.xlabel('Volume $\\AA^3$')
            
# plt.plot(toplotx, to_plot, 'x-', label="RELpsp (Paper lists: 6.432 $\\AA$)")
    
# # Setting the fit label
# fit_label = "Fit (" + str(round(a0_optimized,3)) + " $\\AA$) B$_0$: " + str(round(B,3)) + " GPa, $\\Omega_0$:" + str(round(v0_optimized,3))+ " $\\AA^3$"

# plt.ylabel('Total Energy (eV)')
# plt.plot(fitx, fity, '--', label=fit_label)
# plt.title('Lattice Constant Optimization')
# plt.legend(loc='upper left')
# plt.savefig("./"+dirname+'/LatticeConstant.pdf')

# endTime = time.time()	
# elapsedTimeMin = (endTime - startTime)/60.00
# print("Elapsed Time: " + str(elapsedTimeMin) + " mins")
# file_toWrite_output = open("./"+dirname+"/energies.my", "a")
# file_toWrite_output.write("Elapsed Time: " + str(elapsedTimeMin) + " mins\n")

if __name__ == "__main__":
    run = optimize_parameter()
