"""
This file is intended to replace all run creator files
"""

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
        """

        # Here all the args are set

        # Here all the kwargs are set
        self.run_name = kwargs.get("name", "run")
        self.calculator = kwargs.get("calculator", "espresso")
        self.a0 = kwargs.get("a0", [6.6])
        self.KE_cut = kwargs.get("KE_cut", [20])
        self.k = kwargs.get("k", [2])
        self.pseudopotentials = kwargs.get("pseudopotentials", {'Sn': 'Sn_ONCV_PBE_FR-1.1.upf'})

        # Here goes the PBS init stuff
        self.walltime_mins = 30
        self.nodes = 2
        self.procs = 8

        # Initializations
        self.d = f"^"  # Here you can set the desired delimiter
        self.equals = f"="  # Here you can set the desired symbol for value assigner

        # Some other calc that needs to be done
        

        def make_runs(*args, **kwargs):
            """This method makes all the runs that are required for a single calculation"""
        for KE_cut_i in self.KE_cut:
            for a0_i in self.a0:
                for k_i in self.k:
                    for R_i in self.R:
                        run_name = f"QE{d}KE{KE_cut_i}{d}K{k_i}{d}R{R_i}{d}a{a0_i}{d}PP={self.PP}{d}calc={calc}{d}{}"
                        self.structure.set_cell([(b, 0, 0), (0, b, 0), (0, 0, b)], scale_atoms=True)
                        ase.io.write(f"{run_name}.in", bulk, format = "espresso-in", 

