from ase.io import read, write
from ase.atom import Atom
from ebk import get_machine_paths
import numpy as np

xyz_path = get_machine_paths()["xyz"]
EDT12_path = f"{xyz_path}/1,2-ethaneDithiol_relaxed.xyz"
BDT12_path = f"{xyz_path}/1,2-benzeneDithiol_relaxed.xyz"
BDT14_path = f"{xyz_path}/1,4-benzeneDithiol_relaxed_KE50.xyz"

class Insert_ligand():
    """
    This class is made to be strictly for ligands so that we dont have to go about attaching ligands manually. 
    This class will isert a ligand at a given position and at a given orientation.
    """
    def __init__(self, *args, **kwargs):
        self.name = kwargs.get("name", "EDT12")
        self.path = kwargs.get("path", f"{xyz_path}/Ligands/Relaxed/{self.name}_LDA.xyz")
        self.atoms = read(self.path, index=None, format="xyz")

    def orient(self, direction):
        """
        Plan to make this into a function that takes in  direction and orients the ligands in that direction.
        This might have to be a ligands specific fucntion though.
        """
        if direction == [0,0,1]:
            self.atoms.rotate(90, 'y')
        if direction == [0,0,-1]:
            self.atoms.rotate(-90, 'y')
        if direction == [1,1,0]:
            self.atoms.rotate(45, 'z')
        # if direction == [0,1,0]: # Under construction
        #     self.atoms.rotate(90, 'z')

    def align_atoms_along_axis(self, atom1, atom2, pivot = [0,0,0], axis = [0,0,1]):
        """
        Aligns atom1 and atom 2 to be on the same line (eg: have the same z coordinate)
        What we want to do is to rotate the ligand so that the two atoms are in the same z coordinate
        Inputs:
            atom1: the atom number in the atoms object
            atom2: the atom number in the atoms object
            axis : the axis where we want to align on
            pivot: coordinates of the pivot point.
        """
        atom2 = list(self.atoms[atom2].position)
        atom1 = list(self.atoms[atom1].position)
        # print(f"atom1: {atom1}")
        # print(f"atom2: {atom2}")
        if axis == [0,0,1]:
            y = atom2[1]-atom1[1]
            z = atom2[2]-atom1[2]
            theta_radians = np.arctan(y/z)
            theta = theta_radians*360/(2*np.pi)
            self.atoms.rotate(theta, 'x')

    def make_centre(self, atom1, atom2):
        """Centres the ligand at the bisector of the two atoms."""
        atom2 = self.atoms[atom2].position
        atom1 = self.atoms[atom1].position
        diff_vec = [0, 0, 0]
        for i in range(3):
            diff_vec[i] = (atom1[i] + atom2[i])/2
        
        for atom in self.atoms:
            for i in range(3):
                atom.position[i] = atom.position[i] - diff_vec[i]

    def attach_to_slab(self, slab, sb, sp, lb, lp, center_atom1, center_atom2, bond_length = 0):
        """
        Attaches the ligand by removing the passivation atoms on the slab and ligand.
        Will align the passivation atom vectors of the slab and ligand. There by preserving the dihedral angle.
        Then the slab passivation atom is removed and replced by the ligand non passivation atom at site of attachment.
        The ligands bonding atoms will also be replaced.
        Requirements:
            Slab will need to be passivated
        Inputs,
            slab: will be a zincblende type object which will attach onto the ligand
            sb  : slab bulk atom at site
            sp  : slab passivation atom at site
            lb  : ligand bulk atom at site
            lp  : ligand passivation atom at site
            center_atom1,2: position of the centre atoms of the ligand so that we can set teh 0,0,0 for the system. 
            bond_length = Here we can set the coordinates for the bond length of the attaching bond (eg Sn-S bond). if length is zero no adjustment is made
                This is done by repositioning the slab passivation atom so that it will be at the right distance when being replaced.
        """
        self.lb = lb  # Saving the ligand bulk site globaly
        if bond_length != 0 :
            pass
        Vsp = [0, 0, 0]  # The vector between the slab atoms and the passivation atom
        Vlp = [0, 0, 0]  # The vector between the ligand atoms and the attaching atom
        for i in range(3):
            # Getting the vectors
            Vsp[i] = slab[sb].position[i] - slab[sp].position[i]
            Vlp[i] = self.atoms[lp].position[i] - self.atoms[lb].position[i]
        # self.edit()
        self.atoms.rotate(Vlp, Vsp, center = self.atoms[lb].position)
        diff_vec = [0, 0, 0]
        for i in range(3):
            diff_vec[i] = slab[sp].position[i] - self.atoms[lb].position[i]
        for i,atom in enumerate(slab):
            atom.position = (atom.position[0] - diff_vec[0], atom.position[1] - diff_vec[1], atom.position[2] - diff_vec[2])
        # Deleting the relevant atoms
        del slab[sp]
        del self.atoms[lp]
        for i,atom in enumerate(slab):
            self.atoms.append(atom)
            if i == sb:
                # print(sb)
                self.sb = len(self.atoms)-1  # saving the slab bulk site globally
                # print(self.sb)
        self.make_centre(center_atom1, center_atom2)
        # self.edit()

    def edit(self):
        self.atoms.edit()

    def get_chemical_symbols(self):
        return self.atoms.get_chemical_symbols()

    def get_positions(self):
        return self.atoms.get_positions()

    def list_all_atoms(self):
        syms = self.get_chemical_symbols()
        coor = self.get_positions()
        return [f"{syms[n]}: {coor[n]}" for n in range(len(syms))]

    def save(self, format):
        write(f"{self.name}.{format}", self.atoms)

    def invert_legacy(self):
        for atom in self.atoms:
            for i in range(3):
                atom.position[i] = -atom.position[i]

    def find_cell(self, a0, **kwargs):
        """x, y, z will add extra vacuum"""
        x = kwargs.get("x", 0)  # Any vaccuum that you want to insert
        y = kwargs.get("y", 0)
        z = kwargs.get("z", 0)
        max_x = 0
        max_y = 0
        max_z = 0
        min_x = 1000000
        min_y = 1000000
        min_z = 1000000
        for atom in self.atoms:
            if atom.position[0] > max_x: max_x = atom.position[0]
            if atom.position[1] > max_y: max_y = atom.position[1]
            if atom.position[2] > max_z: max_z = atom.position[2]
            if atom.position[0] < min_x: min_x = atom.position[0]
            if atom.position[1] < min_y: min_y = atom.position[1]
            if atom.position[2] < min_z: min_z = atom.position[2]

        # calculating the cell offsets in teh x, y, z directions. Here we assume the ligand is centred at the origin to startwith
        x_offset = a0/4-self.atoms[self.sb].position[0]*2
        y_offset = -a0/4-self.atoms[self.sb].position[1]*2

        # These are the cell vectors
        x_cell = (a0 + x, 0, 0)
        y_cell = (0, a0 + y, 0)
        z_cell = (x_offset, y_offset, max_z-min_z + a0/4)  # Here we have to add a0/4 since there is no information on that.

        cell = [x_cell, y_cell, z_cell]
        translation_vec = [(x_cell[0]+x_offset)/2, (y_cell[1]+y_offset)/2, (z_cell[2])/2]
        for atom in self.atoms:
            atom.position = (atom.position[0] + translation_vec[0], atom.position[1] + translation_vec[1], atom.position[2] + translation_vec[2])
        self.atoms.set_cell(cell)
        return cell

    def attach(self, attach_to, attach_at, attach_through, atoms_to_delete = [], move_ligand = False):
        """
        This is the main function that attaches the ligand to an atoms type object
        attach_to: The atoms type object that the ligand will attach to
        attach_at: Will attach at this atom in the atoms object(attach_to) and if there is an atom there already it will remove it and attach the ligand
        attach_through: The ligand will remove this atom and will attach through the resulting dangling bond. to the attach_at site in the attach_to object.
        """
        # print(attach_to.atoms)
        # site_attach_at = attach_to.atoms[attach_at]  # Use this if the object is atoms type.
        site_attach_at = attach_to[attach_at]  # Use this if the object is not of atoms type. eg: slabs are bulk type and not atoms type.
        site_attach_through = self.atoms[attach_through]

        diff_vec = [0, 0, 0]
        for i in range(3):
            diff_vec[i] = site_attach_at.position[i] - site_attach_through.position[i]

        # debugging purposes to view atoms
        # attach_to.edit()
        # self.atoms.edit()

        for atom in self.atoms:
            atom.position = (atom.position[0] + diff_vec[0], atom.position[1] + diff_vec[1], atom.position[2] + diff_vec[2])

        #deleting the relevant atoms
        # del attach_to.atoms[attach_at]
        del attach_to[attach_at]
        # write("attach_to.ion.xyz", attach_to.atoms)
        write("attach_to.ion.xyz", attach_to)
        # del self.atoms[attach_through]  # Legacy code where we usually had the atom attached_through being deleted.
        for x in atoms_to_delete:
            del self.atoms[x]
        write(f"{self.name}.ion.xyz", self.atoms)

        # # Logging all the indeces in case we have to constrain them
        # attach_to.NP_atoms = []
        # for atom in attach_to.atoms:
        #     attach_to.NP_atoms.append(atom.index)

        # attach_to.ligand_atoms = []
        # for atom in self.atoms:
        #     # attach_to.atoms.append(atom)
        #     attach_to.append(atom)
        #     attach_to.ligand_atoms.append(len(attach_to)-1)
        #     # print(attach_to.ligand_atoms)
        # return attach_to

        for atom in attach_to:
            self.atoms.append(atom)

        # for atom in attach_to:

    def adsorb(self, attach_to, attach_at, attach_through, z_space = 1, slab_atoms_to_delete = [], ligand_atoms_to_delete = [], move_ligand = False):
        """
        This method in contrast to attach will make the ligand attach straight on by aligning atoms so that it will just sit on the surface. There is also the option to 
        get rid of passivation and end atoms.
        Right now works only in the Z axis. The centre of the ligand is taken to be [0,0,0] and therefore the attach_on atom will have z coordinates of 0
        """
        centre_of_ligand = [0,0,0]  # This is not hard coded and should be changed. Maybe a self.centre_of_ligand property to be added.
        site_attach_at = attach_to[attach_at]  # Use this if the object is not of atoms type. eg: slabs are bulk type and not atoms type.
        site_attach_through = self.atoms[attach_through]  # This is for the ligand

        diff_vec = [0, 0, 0]
        for i in range(3):
            diff_vec[i] = site_attach_at.position[i] - centre_of_ligand[i]

        # Realigning the attach at position to be inline with the centre of the ligand
        z_offset = site_attach_at.position[2] - centre_of_ligand[2]
        for atom in attach_to:
            atom.position = (atom.position[0] - diff_vec[0], atom.position[1] - diff_vec[1], atom.position[2] - z_space)

        for x in ligand_atoms_to_delete:
            del self.atoms[x]

        for x in slab_atoms_to_delete:
            del attach_to[x]

        for atom in attach_to:
            self.atoms.append(atom)

    def make_inversion_symmetric(self):
        """
        This will make the ligand inversion symmetric. 
        Funtionality till now includes only inversion symmetry along teh z axis.
            What this means is that all Z>0 atoms will be rewritten so that we have inversion symmetry
        """
        zg0 = []  # Atoms that are at z>0
        ze0 = []  # Atoms that are at z=0
        zl0 = []  # Atoms that are at z<0

        for i,v in enumerate(self.atoms):
            if v.position[2] > 0: zg0.append(i)
            elif v.position[2] == 0: ze0.append(i)
            elif v.position[2] < 0: zl0.append(i)

        # print(f"zg0: {zg0}")
        # print(f"ze0: {ze0}")
        # print(f"zl0: {zl0}")
        # print(len(self.atoms))

        for x in zl0:
            self.atoms.append(self.atoms[x])
            last_atom = len(self.atoms)-1
            for i in range(3):
                self.atoms[last_atom].position[i] = -self.atoms[last_atom].position[i]

        for x in range(len(self.atoms), -1, -1):
            if x in zg0:
                if x < self.sb: self.sb = self.sb - 1
                elif x == self.sb: self.sb = "Has been deleted"
                del self.atoms[x]

class BDT14(Insert_ligand):
    def __init__(self, *args, **kwargs):
        self.name = kwargs.get("name", "BDT14")
        self.atoms = read(BDT14_path, index=None, format="xyz")
        super().__init__()

class BDT12(Insert_ligand):
    def __init__(self, *args, **kwargs):
        self.name = kwargs.get("name", "BDT12")
        self.atoms = read(BDT12_path, index=None, format="xyz")
        super().__init__()

class EDT12(Insert_ligand):
    def __init__(self, *args, **kwargs):
        self.name = kwargs.get("name", "EDT12")
        self.atoms = read(EDT12_path, index=None, format="xyz")
        super().__init__()


#####################################################################################################################################################
# The code below is the working version as at 2021/01/18.
# This code was finally modified to insert ligands to a surface and also new Functionality like inversion was coded in

# from ase.io import read, write
# from ase.atom import Atom
# from ebk import get_machine_paths

# xyz_path = get_machine_paths()["xyz"]
# EDT12_path = f"{xyz_path}/1,2-ethaneDithiol_relaxed.xyz"
# BDT12_path = f"{xyz_path}/1,2-benzeneDithiol_relaxed.xyz"
# BDT14_path = f"{xyz_path}/1,4-benzeneDithiol_relaxed_KE50.xyz"

# class Insert_ligand():
#     """
#     This class is made to be strictly for ligands so that we dont have to go about attaching ligands manually. 
#     This class will isert a ligand at a given position and at a given orientation.
#     """
#     def __init__(self, *args, **kwargs):
#         self.name = kwargs.get("name", "EDT12")
#         self.path = kwargs.get("path", f"{xyz_path}/Ligands/Relaxed/{self.name}_LDA.xyz")
#         self.atoms = read(self.path, index=None, format="xyz")

#     def orient(self, direction):
#         """
#         Plan to make this into a function that takes in  direction and orients the ligands in that direction.
#         This might have to be a ligands specific fucntion though.
#         """
#         if direction == [0,0,1]:
#             self.atoms.rotate(90, 'y')
#         if direction == [0,0,-1]:
#             self.atoms.rotate(-90, 'y')
#         if direction == [1,1,0]:
#             self.atoms.rotate(45, 'z')

#     def edit(self):
#         self.atoms.edit()

#     def get_chemical_symbols(self):
#         return self.atoms.get_chemical_symbols()

#     def get_positions(self):
#         return self.atoms.get_positions()

#     def list_all_atoms(self):
#         syms = self.get_chemical_symbols()
#         coor = self.get_positions()
#         return [f"{syms[n]}: {coor[n]}" for n in range(len(syms))]

#     def save(self, format):
#         write(f"{self.name}.{format}", self.atoms)

#     def invert(self):
#         for atom in self.atoms:
#             for i in range(3):
#                 atom.position[i] = -atom.position[i]

#     def find_cell(self, **kwargs):
#         """X, y, z will add extra vacuum"""
#         x = kwargs.get("x", 0)
#         y = kwargs.get("y", 0)
#         z = kwargs.get("z", 0)
#         max_x = 0
#         max_y = 0
#         max_z = 0
#         min_x = 1000000
#         min_y = 1000000
#         min_z = 1000000
#         for atom in self.atoms:
#             if atom.position[0] > max_x: max_x = atom.position[0]
#             if atom.position[1] > max_y: max_y = atom.position[1]
#             if atom.position[2] > max_z: max_z = atom.position[2]
#             if atom.position[0] < min_x: min_x = atom.position[0]
#             if atom.position[1] < min_y: min_y = atom.position[1]
#             if atom.position[2] < min_z: min_z = atom.position[2]
#         x_vec = max_x-min_x
#         y_vec = max_y-min_y
#         z_vec = max_z-min_z
#         self.atoms.set_cell([x_vec/2+x, y_vec/2+y, z_vec+z])
#         return [x_vec/2, y_vec/2, z_vec]

#     # def set_cell(self):
#     #     vecs = self.find_cell()

#     def attach(self, attach_to, attach_at, attach_through):
#         """
#         This is the main function that attaches the ligand to an atoms type object
#         attach_to: The atoms type object that the ligand will attach to
#         attach_at: Will attach at this atom in the atoms object(attach_to) and if there is an atom there already it will remove it and attach the ligand
#         attach_through: The ligand will remove this atom and will attach through the resulting dangling bond. to the attach_at site in the attach_to object.
#         """
#         # print(attach_to.atoms)
#         # site_attach_at = attach_to.atoms[attach_at]
#         site_attach_at = attach_to[attach_at]
#         site_attach_through = self.atoms[attach_through]

#         diff_vec = [0, 0, 0]
#         for i in range(3):
#             diff_vec[i] = site_attach_at.position[i] - site_attach_through.position[i]

#         # debugging purposes to view atoms
#         # attach_to.edit()
#         # self.atoms.edit()

#         for atom in self.atoms:
#             atom.position = (atom.position[0] + diff_vec[0], atom.position[1] + diff_vec[1], atom.position[2] + diff_vec[2])

#         #deleting the relevant atoms
#         # del attach_to.atoms[attach_at]
#         del attach_to[attach_at]
#         # write("attach_to.ion.xyz", attach_to.atoms)
#         write("attach_to.ion.xyz", attach_to)
#         del self.atoms[attach_through]
#         write(f"{self.name}.ion.xyz", self.atoms)

#         # # Logging all the indeces in case we have to constrain them
#         # attach_to.NP_atoms = []
#         # for atom in attach_to.atoms:
#         #     attach_to.NP_atoms.append(atom.index)

#         # attach_to.ligand_atoms = []
#         # for atom in self.atoms:
#         #     # attach_to.atoms.append(atom)
#         #     attach_to.append(atom)
#         #     attach_to.ligand_atoms.append(len(attach_to)-1)
#         #     # print(attach_to.ligand_atoms)
#         # return attach_to

#         for atom in attach_to:
#             self.atoms.append(atom)

#         # for atom in attach_to:

# class BDT14(Insert_ligand):
#     def __init__(self, *args, **kwargs):
#         self.name = kwargs.get("name", "BDT14")
#         self.atoms = read(BDT14_path, index=None, format="xyz")
#         super().__init__()

# class BDT12(Insert_ligand):
#     def __init__(self, *args, **kwargs):
#         self.name = kwargs.get("name", "BDT12")
#         self.atoms = read(BDT12_path, index=None, format="xyz")
#         super().__init__()

# class EDT12(Insert_ligand):
#     def __init__(self, *args, **kwargs):
#         self.name = kwargs.get("name", "EDT12")
#         self.atoms = read(EDT12_path, index=None, format="xyz")
#         super().__init__()