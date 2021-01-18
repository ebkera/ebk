from ase.io import read, write
from ase.atom import Atom
from ebk import get_machine_paths

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

    def invert(self):
        for atom in self.atoms:
            for i in range(3):
                atom.position[i] = -atom.position[i]

    def find_cell(self, **kwargs):
        """X, y, z will add extra vacuum"""
        x = kwargs.get("x", 0)
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
        x_vec = max_x-min_x
        y_vec = max_y-min_y
        z_vec = max_z-min_z
        self.atoms.set_cell([x_vec/2+x, y_vec/2+y, z_vec+z])
        return [x_vec/2, y_vec/2, z_vec]

    # def set_cell(self):
    #     vecs = self.find_cell()

    def attach(self, attach_to, attach_at, attach_through):
        """
        This is the main function that attaches the ligand to an atoms type object
        attach_to: The atoms type object that the ligand will attach to
        attach_at: Will attach at this atom in the atoms object(attach_to) and if there is an atom there already it will remove it and attach the ligand
        attach_through: The ligand will remove this atom and will attach through the resulting dangling bond. to the attach_at site in the attach_to object.
        """
        # print(attach_to.atoms)
        # site_attach_at = attach_to.atoms[attach_at]
        site_attach_at = attach_to[attach_at]
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
        del self.atoms[attach_through]
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