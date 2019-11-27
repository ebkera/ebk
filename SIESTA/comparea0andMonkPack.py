import matplotlib
matplotlib.use('Agg')  # no UI backend required if working in the wsl without a UI
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np

# Different x values so that we dont have to do all the points for all the different calculations
a_0 = [6.20,6.22,6.24,6.26,6.28,6.30,6.32,6.34,6.36,6.38,6.40,6.432,6.457,6.475,6.50,6.55,6.60,6.672]
# The data for different calculations
E_1 = [-192.950371,-192.979012,-193.005200,-193.029386,-193.051297,-193.070882,-193.088247,-193.103673,-193.116971,-193.128332,-193.137867,-193.149471,-193.155701,-193.158812,-193.160842,-193.157927,-193.147302,-193.119352]
E_2 = [-195.629572,-195.645151,-195.658392,-195.669392,-195.678244,-195.684983,-195.689658,-195.692276,-195.692809,-195.691379,-195.688262,-195.679847,-195.670365,-195.662027,-195.648466,-195.614450,-195.572594,-195.499551]
E_3 = [-195.804797,-195.820201,-195.833257,-195.844059,-195.852708,-195.859237,-195.863686,-195.866063,-195.866370,-195.864708,-195.861317,-195.852459,-195.842626,-195.834010,-195.820063,-195.785192,-195.742416,-195.667764]
E_4 = [-195.826808,-195.842368,-195.855581,-195.866536,-195.875342,-195.882026,-195.886631,-195.889162,-195.889626,-195.888122,-195.884885,-195.876276,-195.866640,-195.858161,-195.844405,-195.809906,-195.767495,-195.693356]
E_5 = [-195.830961,-195.846604,-195.859900,-195.870940,-195.879831,-195.886602,-195.891294,-195.893914,-195.894467,-195.893053,-195.889906,-195.881444,-195.871922,-195.863524,-195.849878,-195.815583,-195.773338,-195.699403]
E_6 = [-195.831986,-195.847662,-195.860993,-195.872068,-195.880995,-195.887804,-195.892534,-195.895193,-195.895786,-195.894413,-195.891309,-195.882916,-195.873450,-195.865094,-195.851508,-195.817340,-195.775234,-195.701507]
E_7 = [-195.832281,-195.847971,-195.861316,-195.872405,-195.881348,-195.888172,-195.892918,-195.895594,-195.896205,-195.894850,-195.891764,-195.883401,-195.873960,-195.865621,-195.852060,-195.817941,-195.775904,-195.702447]
E_8 = [-195.832375,-195.848071,-195.861421,-195.872516,-195.881465,-195.888296,-195.893049,-195.895732,-195.896351,-195.895004,-195.891927,-195.883581,-195.874153,-195.865825,-195.852279,-195.818194,-195.776192,-195.702756]
E_9 = [-195.832407,-195.848104,-195.861457,-195.872555,-195.881506,-195.888340,-195.893096,-195.895782,-195.896405,-195.895062,-195.891990,-195.883651,-195.874229,-195.865906,-195.852368,-195.818304,-195.776340,-195.702977]
E_10 = [-195.832418,-195.848116,-195.861470,-195.872569,-195.881521,-195.888356,-195.893114,-195.895802,-195.896426,-195.895086,-195.892015,-195.883680,-195.874263,-195.865943,-195.852409,-195.818356,-195.776402,-195.703066]

# THese are the a0,B0, Omega0 values for the energies
E1 = [6.509,39.381,68.917,1]
E2 = [6.357,52.777,64.324,2]
E3 = [6.353,53.467,64.176,3]
E4 = [6.357,53.499,64.176,4]
E5 = [6.357,53.038,64.324,5]
E6 = [6.357,53.016,64.324,6]
E7 = [6.357,53.018,64.324,7]
E8 = [6.357,53.002,64.324,8]
E9 = [6.357,52.994,64.324,9]
E10 = [6.357,52.989,64.324,10]

Energies = []
font = {'size'   : 6}

matplotlib.rc('font', **font)
# The Individual plots.
#plt.plot(a_0, E_1, 'x-', label="a$_0$:" + str(E1[0]) + " $\\AA$ B$_0$: " + str(E1[1]) + " GPa, $\\Omega_0$:" + str(E1[2])+ " $\\AA^3$, MP: "+ str(E1[3]))
#plt.plot(a_0, E_2, 'x-', label="a$_0$:" + str(E2[0]) + " $\\AA$ B$_0$: " + str(E2[1]) + " GPa, $\\Omega_0$:" + str(E2[2])+ " $\\AA^3$, MP: "+ str(E2[3]))
plt.plot(a_0, E_3, 'x-', label="a$_0$:" + str(E3[0]) + " $\\AA$ B$_0$: " + str(E3[1]) + " GPa, $\\Omega_0$:" + str(E3[2])+ " $\\AA^3$, MP: "+ str(E3[3]))
plt.plot(a_0, E_4, 'x-', label="a$_0$:" + str(E4[0]) + " $\\AA$ B$_0$: " + str(E4[1]) + " GPa, $\\Omega_0$:" + str(E4[2])+ " $\\AA^3$, MP: "+ str(E4[3]))
plt.plot(a_0, E_5, 'x-', label="a$_0$:" + str(E5[0]) + " $\\AA$ B$_0$: " + str(E5[1]) + " GPa, $\\Omega_0$:" + str(E5[2])+ " $\\AA^3$, MP: "+ str(E5[3]))
plt.plot(a_0, E_6, 'x-', label="a$_0$:" + str(E6[0]) + " $\\AA$ B$_0$: " + str(E6[1]) + " GPa, $\\Omega_0$:" + str(E6[2])+ " $\\AA^3$, MP: "+ str(E6[3]))
plt.plot(a_0, E_7, 'x-', label="a$_0$:" + str(E7[0]) + " $\\AA$ B$_0$: " + str(E7[1]) + " GPa, $\\Omega_0$:" + str(E7[2])+ " $\\AA^3$, MP: "+ str(E7[3]))
plt.plot(a_0, E_8, 'x-', label="a$_0$:" + str(E8[0]) + " $\\AA$ B$_0$: " + str(E8[1]) + " GPa, $\\Omega_0$:" + str(E8[2])+ " $\\AA^3$, MP: "+ str(E8[3]))
plt.plot(a_0, E_9, 'x-', label="a$_0$:" + str(E9[0]) + " $\\AA$ B$_0$: " + str(E9[1]) + " GPa, $\\Omega_0$:" + str(E9[2])+ " $\\AA^3$, MP: "+ str(E9[3]))
plt.plot(a_0, E_10, 'x-', label="a$_0$:" + str(E10[0]) + " $\\AA$ B$_0$: " + str(E10[1]) + " GPa, $\\Omega_0$:" + str(E10[2])+ " $\\AA^3$, MP: "+ str(E10[3]))

# Setting the fit label
#fit_label = "Fit (" + str(round(a0_optimized,3)) + " $\\AA$) B$_0$: " + str(round(B,3)) + " GPa, $\\Omega_0$:" + str(round(v0_optimized,3))+ " $\\AA^3$"

#plt.plot(fitx, fity, '--', label=fit_label)
plt.ylabel('Total Energy (eV)')
plt.title('Lattice Constant Optimization (Paper: 6.432 $\\AA$)')
plt.xlabel('Lattice Constant a$_0$ ($\\AA$)')
plt.legend(loc='upper left')
plt.savefig('LatticeConstantvsMonk.pdf')