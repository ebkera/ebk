"""
This file is intended to replace all run creator files
"""

from ase.atoms import Atoms

class Calculation_set():
    """
    This class creates an object that will create runs. It will depend entirely on what machine you will run and all that.
    """

    def __init__(self, *args,**kwargs):
        """
        Defined here are all the parameters that you will have to initialize.
        kwargs:
        calculator: string
            Sets the calculator to be used to create runs. Multiple types cannot be set (will have to create two Calculation_set objects).
            Recognized inputs are : espresso, siesta
        structure: [string, [list]]
            sets the type of structure
            Inputs: for string simple_cubic, fcc
                    for 
        structure_parameters: list
            sets the type of structure_parameters. The length of the list depends on type of structure parameter
            Inputs: []
        """

        # Here all the args are set

        # Here all the kwargs are set
        self.description = kwargs.get("description", "run")
        self.calculator = kwargs.get("calculator", "espresso")
        self.a0 = kwargs.get("a0", [6.6])
        self.KE_cut = kwargs.get("KE_cut", [20])
        self.k = kwargs.get("k", [2])
        self.pseudopotentials = kwargs.get("pseudopotentials", {'Sn': 'Sn_ONCV_PBE_FR-1.1.upf'})
        self.R = kwargs.get("R", None)
        self.calc = kwargs.get("calc", "scf")
        self.structure = kwargs.get("structure", ["fcc", [10]])
        self.structure_parameters = [10]

        # Here goes the PBS init stuff
        self.walltime_mins = 30
        self.nodes = 2
        self.procs = 8

        # Initializations
        self.d = f"^"  # Here you can set the desired delimiter
        self.equals = f"="  # Here you can set the desired symbol for value assigner

        # Some other calc that needs to be done
        self.PP = [v+v for k,v in self.pseudopotentials][-1]

    def make_runs(self, *args, **kwargs):
        """This method makes all the runs that are required for a single calculation"""
        for KE_cut_i in self.KE_cut:
            for a0_i in self.a0:
                for k_i in self.k:
                    if self.R == None:
                        self.run_name = f"{self.description}KE{self.equals}{KE_cut_i}{self.d}K{self.equals}{k_i}{self.d}a{self.equals}{a0_i}{self.d}PP{self.equals}{self.PP}{self.d}type{self.equals}{self.calc}"

                        # Creating the run
                        

                        self.atoms_object = Atoms('Sn2', [(0, 0, 0), (0.25, 0.25, 0.25)],  pbc=True)
                        b = a0_i/2.0
                        if self.structure[1].len() == 1:
                            self.atoms_object.set_cell([(0, b, b), (b, 0, b), (b, b, 0)], scale_atoms=True)
                        else:
                            Print("Error in how structure is defined")
                        ase.io.write(f"/{self.run_name}/{run_name}.in", self.atoms_object, format = "espresso-in", 
                                        label           = f"{run_name}",
                                        pseudopotentials= self.pseudopotentials,
                                        pseudo_dir      = "/mnt/c/Users/Eranjan/Desktop/PseudopotentialDatabase",
                                        kpts            = (k_i, k_i, k_i),
                                        ecutwfc         = KE_cut_i,
                                        calculation     = f"{calc}",
                                        lspinorb        = True,
                                        noncolin        = True,
                                        # ecutrho         = KE_cut_i*4,
                                        occupations     = 'smearing',
                                        smearing        = 'gaussian',
                                        degauss         = 0.01,
                                        mixing_beta     = 0.7)
                        # Setting the cell
                        atom_obj.set_cell(cell=[10, 10, 10])

                        # ethaneDithiol = read("../Run_files/XYZdatabase/1,2-benzeneDithiol.xyz", index=None, format="xyz")
                        # ethaneDithiol.set_cell(cell=[10, 10, 10])
                        # self.structure.set_cell([(b, 0, 0), (0, b, 0), (0, 0, b)], scale_atoms=True)
                        # ase.io.write(f"{run_name}.in", bulk, format = "espresso-in", 

if __name__ == "__main__":
    """This is used as an example as to how we use this file."""
    pass    