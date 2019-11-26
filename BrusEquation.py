import os
import matplotlib
matplotlib.use('Agg')  # no UI backend required if working in the wsl without a UI
import sys
import subprocess

import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np

from ebk import *

# Define constants here
Egap = 0.0
epsilonr = 1

def characteristic_reduced_mass(E, r):
    '''E should be in eV'''
    return h**2/(E * 1.6022e-19 * 8 * (r*(1E-9)) **2 )

def brus_print(radius , me, mh):
    '''Prints the Gap energy and the wavelength for a given radius of the dot'''
    Energy = brus(radius, me, mh)
    print(f"Energy gap in eV: {Energy}")
    print(f"Wave length in micro meters: {str(eV2mum(Energy))}")
    # return Energy

def brus(radius, me, mh):
    '''
    Calculates the Gap energy up to the 3rd term of the Brus eq for a given radius of a dot
    | Inputs: radius (float) -  in nm 
    | Outputs energy in eVs
    '''
    Energy = J2eV(Egap + (h**2/(8*(radius)**2)) * (1/me + 1/mh) + 1.786*e**2/(4*np.pi*epsilon0*epsilonr*(radius)*2))
    return Energy

def brus2term_print(radius, me, mh):
    '''Prints the Gap energy and the wavelength for a given radius of the dot'''
    Energy = brus2term(radius, me, mh)
    print(f"Energy gap in eV: {Energy}")
    print(f"Wave length in micro meters: {str(eV2mum(Energy))}")

def brus2term(radius, me, mh):
    '''Calculates the Gap energy up to the 2nd term of the Brus eq for a given radius of a dot'''
    return J2eV(Egap + (h**2/(8*(radius)**2)) * (1/me + 1/mh))

def brus_couloumb_print(radius):
    '''Prints the Gap energy and the wavelength for a given radius of the dot'''
    Energy = brus_couloumb(radius)
    print(f"Energy gap in eV: {Energy}")
    print(f"Wave length in micro meters: {str(eV2mum(Energy))}")
    
def brus_couloumb(radius):
    '''Calcualtes the Gap energy for a given radius of the dot'''
    return J2eV(1.786*(e**2)/(4*np.pi*epsilon0*epsilonr*(radius)))

def Brus2ndvs3dTerms():
    r = np.linspace(1,15,100)
    brus3rd = brus_couloumb(r)
    brus2nd = brus2term(r)

def Oneoverrsquaredunitmaker(oneor):
    return oneor*(1E18)

    # Plotting the figure
    plt.figure()
    plt.plot(r, brus2nd, linewidth=0.4, label="2$^{nd}$ term")
    plt.plot(r, brus3rd, linewidth=0.4, label="3$^{rd}$ term")
    plt.legend(loc='upper right')
    plt.xlabel('Dot radius in nm')
    plt.ylabel('Energy (eV)')
    plt.title('2$^{nd}$ vs 3$^{rd}$ terms of Brus Equation Comparison')
    plt.savefig("Brus2ndvs3rdterm.pdf")

    r = np.linspace(7,18,100)
    brus2nd = brus2term(r)
    brus3rd = brus_couloumb(r)
    plt.figure()
    plt.plot(r, brus2nd, linewidth=0.4, label="2$^{nd}$ term")
    plt.plot(r, brus3rd, linewidth=0.4, label="3$^{rd}$ term")
    plt.legend(loc='upper right')
    plt.xlabel('Dot radius in nm')
    plt.ylabel('Energy (eV)')
    plt.title('2$^{nd}$ vs 3$^{rd}$ terms of Brus Equation Comparison')
    plt.savefig("Brus2ndvs3rdterm7to18.pdf")

def Thirdterm_plotter():
    r = np.linspace(1,15,100)
    y = brus_couloumb(r)

    # Plotting the figure
    plt.figure()
    plt.plot(r, y, linewidth=0.4, label="3$^{rd}$ term")
    plt.legend(loc='upper right')
    plt.xlabel('Dot radius in nm')
    plt.ylabel('Energy of the third term (eV)')
    plt.title('3$^{rd}$ term of Brus Equation ($-1.786 e^2/4 \pi \epsilon_0 \epsilon_r r$)')
    plt.savefig("Brus3rdterm.pdf")

    # Plotting a different range
    r = np.linspace(7,18,100)
    y = brus_couloumb(r)

    # Plotting the figure
    plt.figure()
    plt.plot(r, y, linewidth=0.4, label="3$^{rd}$ term")
    plt.legend(loc='upper right')
    plt.xlabel('Dot radius in nm')
    plt.ylabel('Energy of the third term (eV)')
    plt.title('3$^{rd}$ term of Brus Equation ($-1.786 e^2/4 \pi \epsilon_0 \epsilon_r r$)')
    plt.savefig("Brus3rdterm7to18.pdf")

def different_radii():
    plt.figure()
    plt.ylabel('Gap energy (eV)')
    plt.xlabel('Dot radius in nm')
    plt.plot(r, E, linewidth=0.4, label="For Sn (Brus eq)")
    for x in range(0,len(all_E_points)):
        plt.plot(all_r_points[x], all_E_points[x], 'x', label=all_labels[x])
    plt.legend(loc='upper right')
    # plt.tight_layout() # otherwise the right y-label is slightly clipped
    plt.title("Energy band gap for different radii")
    plt.savefig("different_radii.pdf")

    plt.figure()
    plt.ylabel('Gap energy (eV)')
    plt.xlabel('1/r$^2$ (nm)')
    r2_inverse_Brus = [1/x**2 for x in r] # This is for the Brus equation its has more r points so has to be treated seperately
    plt.plot(r2_inverse_Brus, E, linewidth=0.4, label="For Sn (Brus eq)")
    for x in range(0,len(all_E_points)):
        plt.plot(all_r2_inverse[x], all_E_points[x], 'x', label=all_labels[x])
    plt.legend(loc='upper right')
    # plt.tight_layout() # otherwise the right y-label is slightly clipped
    plt.title("Energy band gap for different radii")
    plt.savefig("different_radii2.pdf")

def different_radii_NoBrus():
    plt.figure()
    plt.ylabel('Gap energy (eV)')
    plt.xlabel('Dot radius (nm)')
    for x in range(0,len(all_E_points)):
        plt.plot(all_r_points[x], all_E_points[x], 'x', label=all_labels[x])
    plt.legend(loc='upper right')
    # plt.tight_layout() # otherwise the right y-label is slightly clipped
    plt.title("Energy band gap for different radii")
    plt.savefig("different_radii_NoBrus.pdf")

    plt.figure()
    plt.ylabel('Gap energy (J)')
    plt.xlabel('1/r$^2$ (1/m$^2$)')
    all_r2_inverse2 = [x*(1E18) for x in all_r2_inverse]
    E2 = [eV2J(x) for x in all_E_points]
    for x in range(0,len(all_E_points)):
        plt.plot(all_r2_inverse2[x], E2[x], 'x', label=all_labels[x])

    # Calculate coefficients for the polynomial
    linear = np.polyfit(all_r2_inverse2, E2, 1) #making the fit
    print("The slope: " + str(linear[0]))
    print("The reduced mass: " + str(h**2/(8*linear[0])))

    #Making the ploynomial
    f_liner = np.poly1d(linear)

    # calculate new x's and y's
    x_new = np.linspace(all_r2_inverse2[0], all_r2_inverse2[-1], 100)
    y_new = f_liner(x_new)

    plt.plot(x_new,y_new, label=" Linear Fit")

    plt.legend(loc='upper left')
    # plt.tight_layout() # otherwise the right y-label is slightly clipped
    plt.title("Energy band gap for different radii")
    plt.savefig("different_radii2_NoBrus.pdf")

r = np.linspace(1,6,100)
r2_inverse = [1/x**2 for x in r]
E = (Egap + (h**2/(8*(r*(1E-9))**2)) * (1/me_Sn + 1/mh_Sn ) + 1.786*(e**2)/(4*np.pi*epsilon0*epsilonr*r*(1E-9)))*6.242e+18

# This is what I got for the non passivated atom
E_mine29 = 0.2
E_x29 = 1.3

# This is what I got for the passivated atoms
E_minePassivated29 = 1.8
E_xPassivated29 = 1.286/2.00

E_minePassivatedadjusted29 = 2.797
E_xPassivatedadjusted29 = 1.286/2.00

E_minePassivated47 = 1.65
E_xPassivated47 = 1.5436/2.00

E_minePassivatedadjusted47 = 2.396
E_xPassivatedadjusted47 = 1.5436/2.00

E_minePassivated71 = 1.36
E_xPassivated71 = 1.60800/2.00

E_minePassivatedadjusted71 = 2.000
E_xPassivatedadjusted71 = 1.60800/2.00

E_minePassivated99 = 1.33
E_xPassivated99 = 1.865/2.00

E_minePassivatedadjusted99 = 1.837838
E_xPassivatedadjusted99 = 1.865/2.00

E_minePassivated123 = 1.2671
E_xPassivated123 = 1.929/2.00

E_minePassivated275 = 0.7417
E_xPassivated275 = 2.5728/2.00

E_minePassivated525 = 0.685
E_xPassivated525 = 3.216/2.00

all_E_points = [E_minePassivated29,E_minePassivatedadjusted29,E_minePassivated47,E_minePassivatedadjusted47,E_minePassivated71,E_minePassivatedadjusted71,E_minePassivated99,E_minePassivatedadjusted99,E_minePassivated123,E_minePassivated275,E_minePassivated525]
all_r_points = [E_xPassivated29,E_xPassivatedadjusted29,E_xPassivated47,E_xPassivatedadjusted47,E_xPassivated71,E_xPassivatedadjusted71,E_xPassivated99,E_xPassivatedadjusted99,E_xPassivated123,E_xPassivated275,E_xPassivated525]
all_labels = ["Sn29 (1.29 nm)","Sn29_relaxed (1.29 nm)","Sn47 (1.54 nm)","Sn47_relaxed (1.54 nm)","Sn71 (1.60 nm)","Sn71_relaxed (1.60 nm)","Sn99 (1.87 nm)","Sn99_relaxed (1.87 nm)","Sn123 (1.93 nm)","Sn275 (2.57 nm)","Sn525 (3.22 nm)"]

# Here we have only the new runs for the optimized passivation bond lengths
# all_E_points = [E_minePassivatedadjusted29,E_minePassivatedadjusted47,E_minePassivatedadjusted71,E_minePassivatedadjusted99]
# all_r_points = [E_xPassivatedadjusted29,E_xPassivatedadjusted47,E_xPassivatedadjusted71,E_xPassivatedadjusted99]
# all_labels = ["Sn29_corrected (1.29 nm)","Sn47_corrected (1.54 nm)","Sn71_corrected (1.60 nm)","Sn99_corrected (1.87 nm)"]

all_r2_inverse = [1/x**2 for x in all_r_points]

# passivated_mass = characteristic_reduced_mass(E_minePassivated29,E_xPassivated29)
# print(passivated_mass)  # (red mass: 5.920557539292997e+31)
# print(characteristic_reduced_mass(11,1.3)) # (red mass: 1.0114285796292201e+31)
# print(1/reduced_mass(me,mh))

# brus_print(8.5)
# brus2term_print(8.5)
# brus_couloumb_print(3)

# for x in range(0,len(all_E_points)):
#     print(characteristic_reduced_mass(all_E_points[x], all_r_points[x]))

# Terms in the 2nd term of the Brus equation
# print(J2eV
#(h**2/(8*(1e-9)**2*reduced_mass(me,mh))))

# Terms in the #rd term of the Brus equation
# print(1.786*(e**2))
# print(1/(4*np.pi*epsilon0*epsilonr))
# print(1.786*(e**2)/(4*np.pi*epsilon0*epsilonr))
# print(1.786*(e**2)/(4*np.pi*epsilon0*epsilonr*(10**(-9))**2))
# print(1.786*(e**2)/(4*np.pi*epsilon0*epsilonr*(10**(-9))))
# brus_couloumb_print(1)

# Energy for the dot band gaps
# print(eV2mum(0.685))

# All the plotters
# Brus2ndvs3dTerms()
# Thirdterm_plotter()
# different_radii()
# different_radii_NoBrus()
# characteristic_reduced_mass(E_Passivated525.685,E_xPassivated525)

# print(f"This is the band diagram {brus_couloumb(1)}")
# brus_print(1)

# print(reduced_mass(1.85340E-32,1.07903E-31))
# print(eV2mum(0.68500))
# different_radii()
# different_radii_NoBrus()