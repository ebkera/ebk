"""
This file is intended to replace all run creator files
"""

class Calculation_set():
    """
    This class creates an object that will create runs. It will depend entirely on what machine you will run and all that.
    """

    def __init__(self, *args,**kwargs):
        """
        Defined here are all the parameters that you will have to set.

        kwargs:
        calculator: string
            Sets the calculator to be used to create runs. Multiple types cannot be set (will have to create two Calculation_set objects).
            Recognized inputs are : esprsso


        """

        # Here all the args are set
        

        # Here all the kwargs are set
        self.calculator = kwargs.get("calculator", "espresso")  

