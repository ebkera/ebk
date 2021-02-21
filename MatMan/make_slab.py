"""
Uses ASE to make slabs:
    ASE does not support higher order surfaces by default but offers the ase.build.cut method. 
    We use this method to quickly make slabs of different surfaces orientations.
"""

# import ase.build.cut as cut
from ase.lattice.compounds import Zincblende  # For future use especially if required for HgTe
from ase.lattice.cubic import Diamond
import ase
import ebk.MatMan

class MakeSlab():
    def __init__(self, *args, **kwargs):
        pass

    def add_vacuum(self, vacuum):
        """ Adds the amount of vacuum in angstroms"""
        ase.build.add_vacuum(self.atoms, vacuum)

    def make_inversion_symmetric(self):
        ebk.MatMan.make_inversion_symmetric(self.atoms)
        self.center_to_cell()

    def center_to_cell(self):
        self.atoms.center()

class Diamond100(MakeSlab):
    def __init__(self, a0, nlayers, *args, **kwargs):
        self.atoms = Diamond(symbol="Sn", latticeconstant=a0, pbc=(1,1,1))
        self.atoms = ase.build.cut(self.atoms, a=(1, 0, 0), b=(0, 1, 0), c=None, clength=None, origo=(0, 0, 0), nlayers=nlayers, extend=1.0, tolerance=0.01, maxatoms=None)
        super().__init__()
