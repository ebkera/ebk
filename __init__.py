"""
This program has all the constants and device parameters that we will need and also some conversion functions that is very commonly used in code in this package
"""

import scipy.constants

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