"""
This file is meant to be for NPs calulations.
"""

from ebk.tb.TB import TB

class ZB(TB):
    """
    This class is for nanoparticles which are ZB like and has an extra atom (like H) for passivation.
    Here naturally we do a single k point calculation.
    """
    def __init__(*args, **kwargs):
        super