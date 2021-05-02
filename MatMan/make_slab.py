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
from ebk.coordinateMaker import CoordinateMaker as i_s_a
from ebk import progress_bar
# from ebk.MatMan.insert_ligands 

class MakeSlab():
    def __init__(self, *args, **kwargs):
        print(kwargs)
        self.a0 = kwargs["a0"]

    def load(self, atoms_object, *args, **kwargs):
        """Loads external structures from atoms like objects"""
        self.atoms = atoms_object

    def add_vacuum(self, vacuum):
        """ Adds the amount of vacuum in angstroms"""
        ase.build.add_vacuum(self.atoms, vacuum)

    def make_inversion_symmetric(self):
        ebk.MatMan.make_inversion_symmetric(self.atoms)
        self.center_to_cell()

    def center_to_cell(self):
        self.atoms.center()

    def passivate_zinc_blende_slab(self, passivation_direction = "z"):
        self.identify_surface_atoms()
        neighbours_list = []
        sign_list = []

        for i in self.surface_atoms_list:
            atom = self.atoms[i]
            # surface_progress.get_progress(i)
            nn_counter = 0
            nearest_neigbours = []
            sign_list_ = []
            for j, atom2 in enumerate(self.atoms):
                if abs(abs(atom.position[0] - atom2.position[0]) - 0.25*self.a0) <= 0.00000000001:
                    if abs(abs(atom.position[1] - atom2.position[1]) - 0.25*self.a0) <= 0.00000000001:
                        if abs(abs(atom.position[2] - atom2.position[2]) - 0.25*self.a0) <= 0.00000000001:
                            nearest_neigbours.append(j)
                            x = atom.position[0] - atom2.position[0]
                            y = atom.position[1] - atom2.position[1]
                            z = atom.position[2] - atom2.position[2]
                            sign_list_.append([x,y,z])
                            nn_counter +=1

                if nn_counter == 4:
                    break

            if nn_counter != 4:
                neighbours_list.append(nearest_neigbours)
                sign_list.append(sign_list_)

        print(neighbours_list)
        print(sign_list)


    def identify_surface_atoms(self):
        """Here we are just recognizing the sites that are not passivatied and have dangling bonds after trimming"""
        print(f"identify_surface_atoms: Started")
        surface_counter = 0
        length = len(self.atoms)
        # surface_progress = progress_bar(length)
        self.surface_atoms_list = []
        for i, atom in enumerate(self.atoms):
            # surface_progress.get_progress(i)
            nearest_neigbours = 0
            for atom2 in self.atoms:
                if abs(abs(atom.position[0] - atom2.position[0]) - 0.25*self.a0) <= 0.00000000001:
                    # print(f"{atom.position[0] - atom2.position[0]} vs {0.25*self.a0}" )
                    if abs(abs(atom.position[1] - atom2.position[1]) - 0.25*self.a0) <= 0.00000000001:
                        if abs(abs(atom.position[2] - atom2.position[2]) - 0.25*self.a0) <= 0.00000000001:
                            nearest_neigbours +=1
                if nearest_neigbours == 4:
                    break

            if nearest_neigbours != 4:
                surface_counter += 1
                self.surface_atoms_list.append(i)

        print(f"identify_surface_atoms: Done: Checked: {i+1} atoms           ")  # White space to prevent ghosting
        print(self.surface_atoms_list)
        return surface_counter

class Diamond_bulk():
    def __init__(self, a0, nlayers, *args, **kwargs):
        self.atoms = Diamond(symbol="Sn", latticeconstant=a0, pbc=(1,1,1))
        super().__init__()

# class Diamond100(MakeSlab):
#     """This works but uses the cut method"""
#     def __init__(self, a0, nlayers, *args, **kwargs):
#         self.atoms = Diamond(symbol="Sn", latticeconstant=a0, pbc=(1,1,1))
#         self.atoms = ase.build.cut(self.atoms, a=(1, 0, 0), b=(0, 1, 0), c=None, clength=None, origo=(0, 0, 0), nlayers=nlayers, extend=1.0, tolerance=0.01, maxatoms=None)
#         super().__init__()

class Diamond100(MakeSlab):
    def __init__(self, a0, nlayers, *args, **kwargs):
        self.atoms = Diamond(symbol="Sn", latticeconstant=a0, directions=[[1,0,0], [0,1,0], [0,0,1]], size=(1,1,nlayers), pbc=(1,1,0), miller=[None, None, [0, 0, 1]])
        super().__init__(**{"a0":a0})

class Diamond110(MakeSlab):
    """Still under construction"""
    def __init__(self, a0, nlayers, *args, **kwargs):
        self.atoms = Diamond(symbol="Sn", latticeconstant=a0, directions=[[1,-1,0], [0,0,-1], [1,1,0]], size=(1,1,nlayers), pbc=(1,1,0), miller=[None, None, [1,1,0]])
        # self.atoms = ase.build.cut(self.atoms, a=(1, 0, 0), b=(0, 1, 0), c=None, clength=None, origo=(0, 0, 0), nlayers=nlayers, extend=1.0, tolerance=0.01, maxatoms=None)
        super().__init__(**{"a0":a0})

class Diamond111(MakeSlab):
    """Creates Slabs in with the 111 direction now lying along the z axis"""
    def __init__(self, a0, nlayers, *args, **kwargs):
        self.atoms = Diamond(symbol="Sn", latticeconstant=a0, directions=[[1,-1,0], [1,1,-2], [1,1,1]], size=(1,1,nlayers), pbc=(1,1,0), miller=[None, None, [1,1,1]])
        super().__init__(**{"a0":a0})

class Diamond210(MakeSlab):
    """Creates Slabs in with the 210 direction now lying along the z axis"""
    def __init__(self, a0, nlayers, *args, **kwargs):
        self.atoms = Diamond(symbol="Sn", latticeconstant=a0, directions=[[1,-2,0], [0,0,-1], [2,1,0]], size=(1,1,nlayers), pbc=(1,1,0), miller=[None, None, [2,1,0]])
        super().__init__(**{"a0":a0})
        
class Diamond211(MakeSlab):
    """Creates Slabs in with the 211 direction now lying along the z axis"""
    def __init__(self, a0, nlayers, *args, **kwargs):
        self.atoms = Diamond(symbol="Sn", latticeconstant=a0, directions=[[1,-2,0], [0,0,-1], [2,1,0]], size=(1,1,nlayers), pbc=(1,1,0), miller=[None, None, [2,1,1]])
        super().__init__(**{"a0":a0})